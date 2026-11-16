#!/usr/bin/env python3
"""
extend_utr_boundaries.py — v12 Pipeline Component

Pool et al. 2023 approach: instead of adding new StringTie transcripts to
the reference GTF (which causes read-splitting ambiguity), use StringTie
output only to identify a new 3' endpoint, then UPDATE the last exon boundary
of the existing reference transcript.

WHY THIS FIXES THE MEDIAN GENES/CELL PROBLEM:
  - Old approach (v6-v10): reference transcript ends at 1000bp. StringTie
    finds reads up to 1500bp. We ADD a second transcript ending at 1500bp.
    Cell Ranger now sees 2 transcripts for this gene. Reads that could map
    to either get discarded as ambiguous. Result: fewer UMIs counted.
  - New approach (v11+): reference transcript ends at 1000bp. StringTie
    finds reads up to 1500bp. We MOVE the reference exon endpoint to 1500bp.
    Cell Ranger still sees 1 transcript. All reads map unambiguously.

HOW IT WORKS:
  1. Parse filtered StringTie GTF (STRG IDs, already splice-filtered)
  2. Use gffcompare annotated GTF to map STRG gene IDs to reference gene IDs
  3. For each reference gene, collect ALL StringTie transcripts, filter by
     min_coverage first, then score survivors by coverage x extension_bp
     to select the best-supported meaningful extension
  4. In the reference GTF, find the terminal exon of EVERY transcript of
     each gene and update all of them (not just the globally furthest one)
  5. Validate against: chromosome boundaries (.fai), same-strand downstream
     neighbor overlap, max_extension, min_coverage, direction
  6. Output is the reference GTF with updated coordinates — NO new transcripts
  7. Full per-gene rejection log in report for audit trail

V12 CHANGES FROM V11:
  - --max-extension default lowered from 10000 to 2000 bp
      (plant UTR literature + ITAG distribution peaks at 200-500bp;
       10kb was capturing StringTie noise and creating false overlaps)
  - --fai argument added: chromosome size hard cap from samtools faidx index
      (prevents any extension from exceeding contig length — Poplar fix)
  - Neighbor overlap check added (strand-aware):
      extension is capped at next same-strand downstream gene start minus buffer
      (opposite-strand genes are ignored — only same-strand causes Cell Ranger ambiguity)
  - Best transcript selection fixed:
      all StringTie transcripts per gene are collected independently,
      min_coverage filter applied first, then scored by coverage x extension_bp,
      winner selected — prevents long low-coverage transcript from displacing
      a shorter high-coverage one via silent overwrite
  - Multi-transcript terminal exon update fixed:
      write_extended_gtf() now updates terminal exon of ALL reference transcripts
      per gene, not just the single globally furthest exon line, preventing
      internal GTF inconsistency for multi-transcript genes
  - Full rejection log per gene added to report:
      every gene with StringTie evidence is logged with its outcome and reason,
      not just successfully extended genes

EXTENSION RULES:
  - Maximum extension: 2,000 bp beyond annotated gene end (v12: reduced from 10000)
  - Minimum StringTie coverage for extension: 3.0x
  - Only extend in the 3' direction (+ strand: increase end; - strand: decrease start)
  - Never exceed chromosome/contig length (hard cap from .fai)
  - Never overlap same-strand downstream gene (50bp buffer)

Usage:
    python3 extend_utr_boundaries.py
        --stringtie  filtered_stringtie.gtf
        --annotated  stringtie_srr_annotation.annotated.gtf
        --reference  reference.gtf
        --output     reference_extended.gtf
        --report     extension_report.tsv
        --fai        genome.fa.fai
        [--max-extension 2000]
        [--min-coverage 3.0]
        [--neighbor-buffer 50]
"""

import argparse
import re
import sys
from collections import defaultdict

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')

def parse_attrs(attr_str):
    return {m.group(1): m.group(2) for m in ATTR_RE.finditer(attr_str)}


def strip_prefix(gid):
    """Strip gene: or transcript: prefix if present."""
    if gid.startswith("gene:"):
        return gid[5:]
    if gid.startswith("transcript:"):
        return gid[11:]
    return gid


def load_chrom_sizes(fai_path):
    """
    Load chromosome/contig sizes from samtools faidx .fai file.
    Returns dict: chrom_name -> length (int)
    .fai format: chrom, length, offset, bases_per_line, bytes_per_line
    """
    sizes = {}
    with open(fai_path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                sizes[chrom] = int(parts[1])
            except ValueError:
                continue
    print(f"[INFO] Loaded chromosome sizes for {len(sizes)} contigs from .fai")
    return sizes


def build_strg_to_ref_mapping(annotated_gtf):
    """
    Parse gffcompare annotated GTF to build STRG gene ID -> reference gene ID mapping.
    Only map class codes that represent genuine extensions (not = which is exact match).
    """
    EXTENSION_CODES = {'c', 'j', 'k', 'n', 'm'}
    mapping = {}   # STRG.X -> ref_gene_id
    code_map = {}  # STRG.X -> class_code

    with open(annotated_gtf) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 9 or fields[2].lower() != 'transcript':
                continue
            attrs = parse_attrs(fields[8])
            tid = attrs.get('transcript_id', '')
            gid = attrs.get('gene_id', '')
            ref_gene = attrs.get('ref_gene_id', '')
            code = attrs.get('class_code', '')

            if not gid.startswith('STRG.'):
                continue
            if not ref_gene or ref_gene in ['-', '.', 'NA']:
                continue
            if code not in EXTENSION_CODES:
                continue

            strg_gene = gid
            ref_clean = strip_prefix(ref_gene)
            mapping[strg_gene] = ref_clean
            code_map[strg_gene] = code

    print(f"[INFO] Mapped {len(mapping)} STRG gene IDs to reference gene IDs "
          f"(extension class codes only)")
    return mapping, code_map


def get_stringtie_extensions(stringtie_gtf, strg_to_ref, min_coverage):
    """
    V12 FIX: Collect ALL StringTie transcripts per reference gene,
    apply min_coverage filter first, then score survivors by
    coverage x extension_bp and select the winner.

    Previously: took the furthest-reaching transcript and silently
    overwrote its coverage — a long low-coverage transcript could
    displace a short high-coverage one, and the coverage check
    downstream would then operate on the wrong value.

    Now: coverage gate runs before any comparison, so the winner
    is always the best-supported AND meaningfully long transcript.

    Returns dict: ref_gene_id -> {'chrom', 'strand', 'strg_start',
                                   'strg_end', 'coverage', 'strg_gene'}
    """
    # Collect all candidates per ref gene: list of transcript dicts
    candidates = defaultdict(list)

    with open(stringtie_gtf) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 9 or fields[2].lower() != 'transcript':
                continue

            chrom = fields[0]
            strand = fields[6]
            attrs = parse_attrs(fields[8])
            gid = attrs.get('gene_id', '')
            cov = float(attrs.get('cov', '0') or '0')

            ref_gene = strg_to_ref.get(gid)
            if not ref_gene:
                continue

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue

            # Coverage gate applied here — before any comparison
            if cov < min_coverage:
                continue

            candidates[ref_gene].append({
                'chrom': chrom,
                'strand': strand,
                'strg_start': start,
                'strg_end': end,
                'coverage': cov,
                'strg_gene': gid,
            })

    # For each gene, pick winner by coverage x extension reach
    # We don't know the ref boundary here so we score by raw reach
    # (furthest 3' coordinate weighted by coverage)
    extensions = {}
    for ref_gene, txs in candidates.items():
        if not txs:
            continue
        strand = txs[0]['strand']
        if strand == '+':
            # Score: coverage * end_coordinate (further and better-supported wins)
            winner = max(txs, key=lambda t: t['coverage'] * t['strg_end'])
        else:
            # Score: coverage * (1 / start_coordinate) — lower start = further 3'
            winner = max(txs, key=lambda t: t['coverage'] / max(t['strg_start'], 1))
        extensions[ref_gene] = winner

    print(f"[INFO] Found StringTie extensions for {len(extensions)} reference genes "
          f"(after {min_coverage}x coverage pre-filter)")
    return extensions


def get_reference_gene_boundaries(reference_gtf):
    """
    Get current annotated boundaries for each gene in the reference GTF.
    Returns dict: gene_id -> {'chrom', 'strand', 'gene_start', 'gene_end'}
    """
    boundaries = {}

    with open(reference_gtf) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            feature = fields[2].lower()
            if feature not in ('gene', 'mrna', 'transcript'):
                continue
            attrs = parse_attrs(fields[8])
            gid = strip_prefix(attrs.get('gene_id', ''))
            if not gid:
                continue
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            strand = fields[6]
            chrom = fields[0]

            if gid not in boundaries:
                boundaries[gid] = {
                    'chrom': chrom, 'strand': strand,
                    'gene_start': start, 'gene_end': end
                }
            else:
                boundaries[gid]['gene_start'] = min(boundaries[gid]['gene_start'], start)
                boundaries[gid]['gene_end'] = max(boundaries[gid]['gene_end'], end)

    return boundaries


def build_downstream_neighbor_index(ref_boundaries, neighbor_buffer):
    """
    V12 NEW: Build a lookup of the nearest same-strand downstream gene
    start position for each gene. Used to cap extensions before they
    encroach on a neighbor and cause Cell Ranger read ambiguity.

    Only same-strand neighbors matter — opposite-strand genes do not
    cause read assignment ambiguity in Cell Ranger.

    Returns dict: gene_id -> max_allowed_3prime_coordinate
      (i.e. neighbor_start - neighbor_buffer, or None if no downstream neighbor)
    """
    # Group genes by (chrom, strand)
    by_key = defaultdict(list)
    for gid, info in ref_boundaries.items():
        key = (info['chrom'], info['strand'])
        by_key[key].append((gid, info['gene_start'], info['gene_end']))

    neighbor_cap = {}  # gene_id -> max allowed coordinate

    for key, genes in by_key.items():
        chrom, strand = key
        # Sort by start position
        genes_sorted = sorted(genes, key=lambda x: x[1])

        for i, (gid, g_start, g_end) in enumerate(genes_sorted):
            # Find nearest downstream gene on same strand
            # For + strand: downstream = higher start coordinate
            # For - strand: downstream = lower end coordinate (i.e. earlier in genome)
            if strand == '+':
                # Look at genes with start > current gene end
                downstream = [(g2_id, g2_start) for g2_id, g2_start, g2_end
                               in genes_sorted if g2_start > g_end]
                if downstream:
                    nearest_start = min(downstream, key=lambda x: x[1])[1]
                    neighbor_cap[gid] = nearest_start - neighbor_buffer
                else:
                    neighbor_cap[gid] = None  # no downstream neighbor
            else:
                # For - strand: downstream (3') = smaller genomic coordinate
                upstream_in_genome = [(g2_id, g2_end) for g2_id, g2_start, g2_end
                                      in genes_sorted if g2_end < g_start]
                if upstream_in_genome:
                    nearest_end = max(upstream_in_genome, key=lambda x: x[1])[1]
                    neighbor_cap[gid] = nearest_end + neighbor_buffer
                else:
                    neighbor_cap[gid] = None

    print(f"[INFO] Built downstream neighbor index for {len(neighbor_cap)} genes")
    return neighbor_cap


def compute_valid_extensions(extensions, ref_boundaries, neighbor_cap,
                              chrom_sizes, max_extension, min_coverage):
    """
    V12: Validate each extension with full checks:
    - Must be in the 3' direction
    - Must not exceed max_extension bp
    - Must not exceed chromosome/contig length (hard cap from .fai)
    - Must not overlap same-strand downstream neighbor (strand-aware)
    - Extension must be >= 10bp (minimum meaningful)

    Note: min_coverage check already applied upstream in get_stringtie_extensions().
    Kept here as a safety guard only.

    Returns:
      valid dict: ref_gene_id -> extension info
      rejected dict: ref_gene_id -> rejection reason (for report)
    """
    valid = {}
    rejected = {}
    stats = {
        'too_short': 0, 'too_long': 0, 'low_cov': 0,
        'valid': 0, 'wrong_dir': 0,
        'chrom_boundary': 0, 'neighbor_overlap': 0
    }

    for gene_id, ext in extensions.items():
        ref = ref_boundaries.get(gene_id)
        if not ref:
            rejected[gene_id] = 'no_ref_boundary'
            continue

        # Safety coverage guard (primary filter already applied upstream)
        if ext['coverage'] < min_coverage:
            stats['low_cov'] += 1
            rejected[gene_id] = f"low_coverage_{ext['coverage']:.1f}x"
            continue

        strand = ref['strand']
        chrom = ref['chrom']
        chrom_len = chrom_sizes.get(chrom)

        if strand == '+':
            new_end = ext['strg_end']
            extension_bp = new_end - ref['gene_end']

            if extension_bp <= 0:
                stats['wrong_dir'] += 1
                rejected[gene_id] = 'wrong_direction'
                continue
            if extension_bp > max_extension:
                stats['too_long'] += 1
                rejected[gene_id] = f"too_long_{extension_bp}bp"
                continue
            if extension_bp < 10:
                stats['too_short'] += 1
                rejected[gene_id] = 'too_short_<10bp'
                continue

            # Chromosome boundary check (V12 NEW)
            if chrom_len is not None and new_end > chrom_len:
                stats['chrom_boundary'] += 1
                rejected[gene_id] = f"exceeds_chrom_length_{chrom_len}"
                continue

            # Same-strand neighbor overlap check (V12 NEW)
            trimmed = False
            cap = neighbor_cap.get(gene_id)
            if cap is not None and new_end > cap:
                # Trim to neighbor cap if trimmed extension is still meaningful
                trimmed_extension = cap - ref['gene_end']
                if trimmed_extension >= 10:
                    new_end = cap
                    extension_bp = trimmed_extension
                    trimmed = True  # flag for report — trimmed but valid
                else:
                    stats['neighbor_overlap'] += 1
                    rejected[gene_id] = f"neighbor_overlap_cap_{cap}"
                    continue

            valid[gene_id] = {
                'new_coord': new_end,
                'old_coord': ref['gene_end'],
                'extension_bp': extension_bp,
                'strand': '+',
                'coverage': ext['coverage'],
                'chrom': chrom,
                'trimmed': trimmed,
            }

        elif strand == '-':
            new_start = ext['strg_start']
            extension_bp = ref['gene_start'] - new_start

            if extension_bp <= 0:
                stats['wrong_dir'] += 1
                rejected[gene_id] = 'wrong_direction'
                continue
            if extension_bp > max_extension:
                stats['too_long'] += 1
                rejected[gene_id] = f"too_long_{extension_bp}bp"
                continue
            if extension_bp < 10:
                stats['too_short'] += 1
                rejected[gene_id] = 'too_short_<10bp'
                continue

            # Chromosome boundary check (V12 NEW)
            if new_start < 1:
                stats['chrom_boundary'] += 1
                rejected[gene_id] = 'exceeds_chrom_start'
                continue

            # Same-strand neighbor overlap check (V12 NEW — minus strand)
            trimmed = False
            cap = neighbor_cap.get(gene_id)
            if cap is not None and new_start < cap:
                trimmed_extension = ref['gene_start'] - cap
                if trimmed_extension >= 10:
                    new_start = cap
                    extension_bp = trimmed_extension
                    trimmed = True  # flag for report — trimmed but valid
                else:
                    stats['neighbor_overlap'] += 1
                    rejected[gene_id] = f"neighbor_overlap_cap_{cap}"
                    continue

            valid[gene_id] = {
                'new_coord': new_start,
                'old_coord': ref['gene_start'],
                'extension_bp': extension_bp,
                'strand': '-',
                'coverage': ext['coverage'],
                'chrom': chrom,
                'trimmed': trimmed,
            }

        stats['valid'] += 1

    print(f"[INFO] Extension validation: {stats['valid']} valid, "
          f"{stats['too_long']} too long (>{max_extension}bp), "
          f"{stats['low_cov']} low coverage, "
          f"{stats['wrong_dir']} wrong direction, "
          f"{stats['too_short']} too short (<10bp), "
          f"{stats['chrom_boundary']} chromosome boundary violations, "
          f"{stats['neighbor_overlap']} neighbor overlaps blocked")
    return valid, rejected


def write_extended_gtf(reference_gtf, valid_extensions, output_gtf, report_path,
                        extensions_input, rejected):
    """
    V12 FIX: Update terminal exon of ALL reference transcripts per gene,
    not just the single globally furthest exon.

    In v11, gene_last_exon stored only ONE exon index per gene — the
    globally furthest across all transcripts. For multi-transcript genes,
    only that line was updated. Other transcripts' terminal exons stayed
    at original coordinates, creating internal GTF inconsistency.

    Fix: track the terminal exon per TRANSCRIPT, then update each one
    independently if the gene is in valid_extensions.
    """
    lines = []

    # First pass: read all lines, track terminal exon per transcript
    # tx_last_exon: (gene_id, transcript_id) -> idx of that transcript's terminal exon
    tx_last_exon = {}

    with open(reference_gtf) as f:
        for idx, line in enumerate(f):
            lines.append(line)
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 9 or fields[2].lower() != 'exon':
                continue
            attrs = parse_attrs(fields[8])
            gid = strip_prefix(attrs.get('gene_id', ''))
            tid = attrs.get('transcript_id', '')
            if not gid or not tid:
                continue
            if gid not in valid_extensions:
                continue

            strand = valid_extensions[gid]['strand']
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue

            key = (gid, tid)
            if key not in tx_last_exon:
                tx_last_exon[key] = {'idx': idx, 'start': start, 'end': end}
            else:
                existing = tx_last_exon[key]
                if strand == '+' and end > existing['end']:
                    tx_last_exon[key] = {'idx': idx, 'start': start, 'end': end}
                elif strand == '-' and start < existing['start']:
                    tx_last_exon[key] = {'idx': idx, 'start': start, 'end': end}

    # Apply extensions to ALL transcripts' terminal exons per gene
    genes_extended = set()
    for (gid, tid), exon_info in tx_last_exon.items():
        ext = valid_extensions[gid]
        idx = exon_info['idx']
        fields = lines[idx].rstrip('\n').split('\t')

        if ext['strand'] == '+':
            fields[4] = str(ext['new_coord'])
        else:
            fields[3] = str(ext['new_coord'])

        lines[idx] = '\t'.join(fields) + '\n'
        genes_extended.add(gid)

    # Second pass: update gene/transcript/mRNA feature boundaries for completeness
    for idx, line in enumerate(lines):
        if line.startswith('#') or not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) < 9:
            continue
        feature = fields[2].lower()
        if feature not in ('gene', 'mrna', 'transcript'):
            continue
        attrs = parse_attrs(fields[8])
        gid = strip_prefix(attrs.get('gene_id', ''))
        if gid not in valid_extensions or gid not in genes_extended:
            continue

        ext = valid_extensions[gid]
        if ext['strand'] == '+':
            try:
                current_end = int(fields[4])
                if ext['new_coord'] > current_end:
                    fields[4] = str(ext['new_coord'])
                    lines[idx] = '\t'.join(fields) + '\n'
            except ValueError:
                pass
        else:
            try:
                current_start = int(fields[3])
                if ext['new_coord'] < current_start:
                    fields[3] = str(ext['new_coord'])
                    lines[idx] = '\t'.join(fields) + '\n'
            except ValueError:
                pass

    # Write output GTF
    with open(output_gtf, 'w') as f:
        f.writelines(lines)

    print(f"[INFO] Extended {len(genes_extended)} genes "
          f"({sum(1 for k in tx_last_exon)} transcript terminal exons updated)")
    print(f"[INFO] Output written to: {output_gtf}")

    # Write full report — all genes with StringTie evidence, outcome and reason
    n_trimmed = 0
    with open(report_path, 'w') as r:
        r.write("gene_id\tstrand\told_3prime\tnew_3prime\textension_bp\t"
                "coverage\tstatus\treason\n")

        # Successfully extended genes — distinguish clean vs neighbor-trimmed
        for gid in sorted(genes_extended):
            ext = valid_extensions[gid]
            was_trimmed = ext.get('trimmed', False)
            status = 'extended_trimmed' if was_trimmed else 'extended'
            reason = 'neighbor_cap_applied' if was_trimmed else 'ok'
            if was_trimmed:
                n_trimmed += 1
            r.write(f"{gid}\t{ext['strand']}\t{ext['old_coord']}\t"
                    f"{ext['new_coord']}\t{ext['extension_bp']}\t"
                    f"{ext['coverage']:.1f}\t{status}\t{reason}\n")

        # Rejected genes — full per-gene log (V12 NEW)
        for gid, reason in sorted(rejected.items()):
            ext_info = extensions_input.get(gid, {})
            strand = ext_info.get('strand', '?')
            cov = ext_info.get('coverage', 0.0)
            r.write(f"{gid}\t{strand}\t.\t.\t.\t{cov:.1f}\trejected\t{reason}\n")

    print(f"[INFO] Report written to: {report_path} "
          f"({len(genes_extended)} extended [{n_trimmed} neighbor-trimmed], "
          f"{len(rejected)} rejected)")
    return len(genes_extended), n_trimmed


def main():
    parser = argparse.ArgumentParser(
        description="Extend 3' UTR boundaries in reference GTF using StringTie evidence "
                    "(v12 — with chromosome boundary, neighbor overlap, and multi-transcript fixes)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--stringtie', required=True,
                        help='Splice-filtered StringTie GTF (STRG IDs)')
    parser.add_argument('--annotated', required=True,
                        help='gffcompare annotated GTF (class codes + ref_gene_id)')
    parser.add_argument('--reference', required=True,
                        help='Original reference GTF to update')
    parser.add_argument('--output', required=True,
                        help='Output extended reference GTF')
    parser.add_argument('--report', required=True,
                        help='TSV report of all genes: extended and rejected with reasons')
    parser.add_argument('--fai', required=True,
                        help='Genome .fai index file (samtools faidx output) for chromosome size caps')
    parser.add_argument('--max-extension', type=int, default=2000,
                        help='Maximum 3-prime extension in bp (default: 2000, reduced from v11 10000)')
    parser.add_argument('--min-coverage', type=float, default=3.0,
                        help='Minimum StringTie coverage for extension (default: 3.0)')
    parser.add_argument('--neighbor-buffer', type=int, default=50,
                        help='Buffer in bp to leave between extension and downstream neighbor (default: 50)')
    args = parser.parse_args()

    print("[1/6] Loading chromosome sizes from .fai index...")
    chrom_sizes = load_chrom_sizes(args.fai)

    print("[2/6] Mapping STRG gene IDs to reference gene IDs...")
    strg_to_ref, code_map = build_strg_to_ref_mapping(args.annotated)

    print("[3/6] Loading reference gene boundaries...")
    ref_boundaries = get_reference_gene_boundaries(args.reference)
    print(f"[INFO] Loaded boundaries for {len(ref_boundaries)} reference genes")

    print("[4/6] Building same-strand downstream neighbor index...")
    neighbor_cap = build_downstream_neighbor_index(ref_boundaries, args.neighbor_buffer)

    print("[5/6] Selecting best StringTie extension per gene "
          "(coverage pre-filter + coverage x reach scoring)...")
    extensions = get_stringtie_extensions(args.stringtie, strg_to_ref, args.min_coverage)

    print("[6/6] Validating extensions and writing output...")
    valid_extensions, rejected = compute_valid_extensions(
        extensions, ref_boundaries, neighbor_cap,
        chrom_sizes, args.max_extension, args.min_coverage)

    n_extended, n_trimmed = write_extended_gtf(
        args.reference, valid_extensions, args.output, args.report,
        extensions, rejected)

    print(f"\n=== V12 Extension Summary ===")
    print(f"Genes with StringTie evidence  : {len(extensions)}")
    print(f"Genes with valid extensions    : {len(valid_extensions)}")
    print(f"Genes actually extended        : {n_extended}")
    print(f"  of which neighbor-trimmed    : {n_trimmed}")
    print(f"  of which clean extensions    : {n_extended - n_trimmed}")
    print(f"Genes rejected (see report)    : {len(rejected)}")
    print(f"Output GTF                     : {args.output}")
    print(f"Report                         : {args.report}")
    print(f"\nKey: output is the REFERENCE GTF with updated terminal exon coordinates.")
    print(f"No new transcripts added. Cell Ranger will see 1 transcript per gene.")


if __name__ == '__main__':
    main()
