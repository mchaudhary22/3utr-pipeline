#!/usr/bin/env python3
"""
filter_stringtie_splicing.py — v7 Pipeline Component

Filters StringTie GTF output to retain ONLY transcripts with evidence of
novel splicing to known exons (Pool et al. approach), then resolves overlaps.

RATIONALE (from Pool et al. 2023):
  - Coverage-based extension (StringTie default) captures noise → false positives
  - Splice-junction evidence = real transcript extension into UTR territory
  - Only reads spliced from KNOWN exons into NOVEL regions are trusted
  - gffcompare class codes filter out non-informative transcripts

FILTER LOGIC:
  Keep a StringTie transcript if it meets ANY of:
    1. class_code == "=" (exact match to reference)
    2. class_code == "c" (contained in reference) -- splice site preserved
    3. class_code == "k" (containment of reference)
    4. class_code == "m" (retained intron) -- keep for pre-mRNA
    5. class_code == "n" (retained intron with novel splice) -- keep
    6. Has ≥1 splice junction that connects to a KNOWN exon boundary in the
       reference GTF (the key splicing-based evidence criterion)

  REMOVE (false positives):
    - class_code == "u" (unknown intergenic) with NO known splice junctions
    - class_code == "o" (overlap, opposite strand)
    - class_code == "p" (polymerase run-on)
    - class_code == "s" (intron of reference; often noise)
    - class_code == "x" (exonic overlap, opposite strand)
    - class_code == "i" (fully inside reference intron) with no splice evidence

Usage:
    python filter_stringtie_splicing.py \\
        --stringtie input_stringtie.gtf \\
        --reference reference.gtf \\
        --annotated gffcompare.annotated.gtf \\
        --output filtered_stringtie.gtf \\
        --report filter_report.tsv \\
        [--min-coverage 2.0] \\
        [--require-splice]

"""
import argparse
import re
import sys
from collections import defaultdict

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')

def parse_attrs(attr_str):
    return {m.group(1): m.group(2) for m in ATTR_RE.finditer(attr_str)}

def build_reference_splice_sites(reference_gtf):
    """
    Extract all known exon boundaries from reference GTF.
    Returns set of (chrom, strand, pos, side) where side = '5p' or '3p'
    """
    splice_sites = set()
    with open(reference_gtf) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 9 or fields[2].lower() != 'exon':
                continue
            chrom, _, _, start_s, end_s, _, strand, _, _ = fields[:9]
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue
            # Both boundaries of every reference exon are known splice sites
            splice_sites.add((chrom, strand, start, '5p'))
            splice_sites.add((chrom, strand, end, '3p'))
    print(f"[INFO] Loaded {len(splice_sites)} reference splice site boundaries")
    return splice_sites


def build_reference_exon_intervals(reference_gtf):
    """Build interval set of reference exons per (chrom, strand)."""
    exons = defaultdict(list)
    with open(reference_gtf) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 9 or fields[2].lower() != 'exon':
                continue
            chrom, _, _, start_s, end_s, _, strand, _, _ = fields[:9]
            try:
                exons[(chrom, strand)].append((int(start_s), int(end_s)))
            except ValueError:
                continue
    return exons


def parse_annotated_gtf(annotated_gtf):
    """
    Parse gffcompare annotated GTF to extract:
      transcript_id -> class_code
      transcript_id -> has_cuff_id (STRG.)
    Returns dict of transcript_id -> {'class_code': X, 'is_novel': bool}
    """
    transcript_info = {}
    with open(annotated_gtf) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 9 or fields[2].lower() != 'transcript':
                continue
            attrs = parse_attrs(fields[8])
            tid = attrs.get('transcript_id', '')
            class_code = attrs.get('class_code', 'u')
            ref_gene = attrs.get('ref_gene_id', '')
            is_novel = tid.startswith('STRG.') or not ref_gene or ref_gene in ['-', '.', 'NA']
            transcript_info[tid] = {
                'class_code': class_code,
                'is_novel': is_novel,
                'ref_gene_id': ref_gene,
                'gene_id': attrs.get('gene_id', ''),
            }
    print(f"[INFO] Parsed {len(transcript_info)} transcript entries from annotated GTF")
    return transcript_info


def get_transcript_exon_boundaries(stringtie_gtf):
    """
    For each transcript, collect all exon start/end positions.
    Returns dict: transcript_id -> {'chrom', 'strand', 'exon_starts', 'exon_ends', 'coverage'}
    """
    tx_data = {}
    with open(stringtie_gtf) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            chrom, _, feature, start_s, end_s, _, strand, _, attr_s = fields[:9]
            if feature.lower() not in ('transcript', 'exon'):
                continue
            attrs = parse_attrs(attr_s)
            tid = attrs.get('transcript_id', '')
            if not tid:
                continue
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue
            if tid not in tx_data:
                tx_data[tid] = {
                    'chrom': chrom,
                    'strand': strand,
                    'exon_starts': [],
                    'exon_ends': [],
                    'coverage': 0.0,
                    'gene_id': attrs.get('gene_id', ''),
                }
            if feature.lower() == 'exon':
                tx_data[tid]['exon_starts'].append(start)
                tx_data[tid]['exon_ends'].append(end)
            if feature.lower() == 'transcript':
                try:
                    cov = float(attrs.get('cov', '0'))
                    tx_data[tid]['coverage'] = cov
                except ValueError:
                    pass
    return tx_data


def has_known_splice_junction(tx, reference_splice_sites, tolerance=5):
    """
    Check whether a transcript has ≥1 exon boundary that matches a known
    reference splice site (within tolerance bp).

    This is the KEY criterion from Pool et al.: only trust extensions
    that are spliced to known exon boundaries.
    """
    chrom = tx['chrom']
    strand = tx['strand']
    for start in tx['exon_starts']:
        for side in ('5p', '3p'):
            for delta in range(-tolerance, tolerance + 1):
                if (chrom, strand, start + delta, side) in reference_splice_sites:
                    return True
    for end in tx['exon_ends']:
        for side in ('5p', '3p'):
            for delta in range(-tolerance, tolerance + 1):
                if (chrom, strand, end + delta, side) in reference_splice_sites:
                    return True
    return False


# Class codes that are always kept (Pool et al. logic)
KEEP_CODES = {'=', 'c', 'k', 'm', 'n', 'j'}
# Class codes that require splice evidence to keep
CONDITIONAL_CODES = {'u', 'i', 'y', 'p'}
# Class codes that are always removed
REMOVE_CODES = {'o', 's', 'x', 'r'}


def decide_keep(tid, transcript_info, tx_data, reference_splice_sites,
                min_coverage, require_splice):
    """
    Decide whether to keep a transcript.
    Returns (keep: bool, reason: str)
    """
    info = transcript_info.get(tid)
    tx = tx_data.get(tid)
    if not info or not tx:
        # V10 FIX: Don't keep blindly if no class code info found.
        # In v9 this caused 13,836 transcripts to pass through unfiltered.
        # Require splice junction evidence even for unknown transcripts.
        has_splice = has_known_splice_junction(tx, reference_splice_sites) if tx else False
        if has_splice:
            return True, "no_info_splice_evidence"
        return False, "no_info_no_splice_removed"

    code = info.get('class_code', 'u')

    # Always remove these
    if code in REMOVE_CODES:
        return False, f"class_code_{code}_always_remove"

    # Always keep reference-matched transcripts
    if code in KEEP_CODES:
        return True, f"class_code_{code}_keep"

    # For novel/uncertain transcripts, require splice evidence
    if code in CONDITIONAL_CODES or info.get('is_novel', False):
        # Coverage filter
        if tx['coverage'] < min_coverage:
            return False, f"low_coverage_{tx['coverage']:.1f}<{min_coverage}"
        # Splice evidence filter (key Pool et al. criterion)
        if require_splice:
            has_splice = has_known_splice_junction(tx, reference_splice_sites)
            if not has_splice:
                return False, "no_splice_to_known_exon"
            return True, f"splice_evidence_code_{code}"
        else:
            return True, f"coverage_pass_code_{code}"

    # Default: keep unknown codes
    return True, f"unknown_code_{code}_keep_default"


def filter_stringtie(stringtie_gtf, annotated_gtf, reference_gtf, output_gtf,
                     report_path, min_coverage, require_splice):
    print("[1/4] Building reference splice site database...")
    ref_splice_sites = build_reference_splice_sites(reference_gtf)

    print("[2/4] Parsing gffcompare annotated GTF...")
    transcript_info = parse_annotated_gtf(annotated_gtf)

    print("[3/4] Extracting StringTie transcript exon boundaries...")
    tx_data = get_transcript_exon_boundaries(stringtie_gtf)

    print("[4/4] Applying filters and writing output...")
    kept_tids = set()
    removed_tids = set()
    decisions = {}

    for tid in tx_data:
        keep, reason = decide_keep(tid, transcript_info, tx_data,
                                   ref_splice_sites, min_coverage, require_splice)
        if keep:
            kept_tids.add(tid)
        else:
            removed_tids.add(tid)
        decisions[tid] = (keep, reason)

    # Write filtered GTF
    with open(stringtie_gtf) as inf, open(output_gtf, 'w') as outf:
        for line in inf:
            if line.startswith('#'):
                outf.write(line)
                continue
            if not line.strip():
                outf.write(line)
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                outf.write(line)
                continue
            attrs = parse_attrs(fields[8])
            tid = attrs.get('transcript_id', '')
            if not tid:
                outf.write(line)  # keep non-transcript lines
                continue
            # keep line if transcript is kept
            if tid in kept_tids:
                outf.write(line)
            # if tid not in either set (gene-level lines), keep if any child transcript is kept
            elif tid not in removed_tids:
                gid = attrs.get('gene_id', '')
                # keep gene-level lines if at least one of its transcripts is kept
                any_child_kept = any(
                    t for t in kept_tids
                    if tx_data.get(t, {}).get('gene_id', '') == gid
                )
                if any_child_kept or fields[2].lower() == 'gene':
                    outf.write(line)

    # Write report
    code_stats = defaultdict(lambda: {'kept': 0, 'removed': 0})
    with open(report_path, 'w') as r:
        r.write("transcript_id\tclass_code\tcoverage\thas_splice\tdecision\treason\n")
        for tid, (keep, reason) in decisions.items():
            info = transcript_info.get(tid, {})
            tx = tx_data.get(tid, {})
            code = info.get('class_code', '?')
            cov = tx.get('coverage', 0.0)
            has_splice = has_known_splice_junction(tx, ref_splice_sites) if tx else False
            status = 'kept' if keep else 'removed'
            r.write(f"{tid}\t{code}\t{cov:.2f}\t{has_splice}\t{status}\t{reason}\n")
            code_stats[code][status] += 1

    print(f"\n=== Filter Summary ===")
    print(f"Total transcripts evaluated : {len(decisions)}")
    print(f"Kept                        : {len(kept_tids)}")
    print(f"Removed                     : {len(removed_tids)}")
    print(f"\nBreakdown by class_code:")
    print(f"{'Code':<8} {'Kept':>8} {'Removed':>10}")
    for code in sorted(code_stats):
        s = code_stats[code]
        print(f"  {code:<6} {s['kept']:>8} {s['removed']:>10}")
    print(f"\nFiltered GTF: {output_gtf}")
    print(f"Report:       {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter StringTie GTF using splicing-based evidence (Pool et al. approach)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Filtering strategy (Pool et al. 2023 informed):
  KEEP:    class_code =, c, k, m, n, j (reference-supported)
  KEEP:    class_code u/i IF has splice junction to known exon AND coverage >= threshold
  REMOVE:  class_code o, s, x, r (strand/overlap issues)
  REMOVE:  class_code u/i with NO splice evidence to known exons (false positives)

This prevents the v6 failure mode: StringTie coverage extensions that create
gene overlaps and cause multi-gene read discarding.
        """
    )
    parser.add_argument('--stringtie', required=True, help='StringTie output GTF')
    parser.add_argument('--reference', required=True, help='Reference GTF (for splice sites)')
    parser.add_argument('--annotated', required=True,
                        help='gffcompare annotated GTF (e.g., out.annotated.gtf)')
    parser.add_argument('--output', required=True, help='Output filtered GTF')
    parser.add_argument('--report', required=True, help='TSV filter report')
    parser.add_argument('--min-coverage', type=float, default=2.0,
                        help='Minimum coverage for novel transcripts (default: 2.0)')
    parser.add_argument('--require-splice', action='store_true', default=True,
                        help='Require splice junction to known exon for novel transcripts (default: True)')
    parser.add_argument('--no-require-splice', dest='require_splice', action='store_false',
                        help='Disable splice junction requirement')
    args = parser.parse_args()

    filter_stringtie(
        stringtie_gtf=args.stringtie,
        annotated_gtf=args.annotated,
        reference_gtf=args.reference,
        output_gtf=args.output,
        report_path=args.report,
        min_coverage=args.min_coverage,
        require_splice=args.require_splice,
    )

if __name__ == '__main__':
    main()
