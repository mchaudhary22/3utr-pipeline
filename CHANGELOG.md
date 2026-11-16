# Changelog

## v1.0.0: first public release (internal development version v13)

First public release of UTRex. It replaces the earlier v1 pipeline code in this repository.

- Pipeline: GFF3 → GTF, same-strand overlap removal, baseline Cell Ranger run, reference-guided StringTie, gffcompare, splice-junction-aware filter, 3′ end extension of the reference genes, final Cell Ranger run.
- The reference GTF keeps one model per gene. UTRex moves the gene's 3′ end; it does not add StringTie transcripts.
- New scripts: `remove_overlap.py`, `filter_stringtie_splicing.py`, `extend_utr_boundaries.py`.
- Removed: `update_ref_gene_format.py` (no longer used).
- Cell Ranger 10.0.0. The Dockerfile now expects a Cell Ranger tarball that you download from 10x Genomics yourself.
- New run options: `-f` (force cell count). `-c` no longer has a default; without `-c` or `-f`, Cell Ranger detects the cell count itself.
- New reports: overlap report, splice filter report, extension report.
- Docker image: `mchaudhary22/3utr-pipeline-v13:latest`.

## Development notes

A short summary of the internal versions that led to v1.0.0. More detail is in the version comments inside the code.

- **v6: overlap removal.** Genes that overlap on the same strand are grouped, and one gene per group is kept (`remove_overlap.py`). Overlapping genes on the same strand make Cell Ranger discard reads as multi-gene.
- **v7: splice-junction filter.** StringTie transcripts are filtered with gffcompare class codes and exon-boundary evidence (`filter_stringtie_splicing.py`), following the approach of Pool et al. 2023.
- **v10: gene-ID fix and filter order.**
  - `gene:` and `transcript:` prefixes are stripped after gffread. Before this, the same gene could appear under two IDs (`gene:X` from gffread and `X` from StringTie), which doubled genes and broke overlap detection.
  - gffcompare and the splice filter now run on the original StringTie (STRG) IDs, before any ID remapping. In v9 the remapping ran first, so many transcripts had no class code and passed the filter unchecked.
  - Transcripts with no class-code information now need exon-boundary evidence to be kept.
- **v11: move the gene end instead of adding transcripts.** Earlier versions added StringTie transcripts to the reference, so one gene could have several overlapping models. Cell Ranger treats reads that fit more than one model as ambiguous, so fewer reads were counted. From v11 on, UTRex moves the end of the existing terminal exon and adds no new transcripts. The final Cell Ranger reference is built from this extended reference GTF.
- **v12: chromosome cap and neighbor check.**
  - `samtools faidx` gives chromosome lengths, and no extension may run past the chromosome end.
  - An extension must stop 50 bp before the next downstream gene on the same strand. It is trimmed if at least 10 bp of extension is left, and rejected otherwise.
  - The coverage filter runs before the best StringTie transcript is chosen for each gene.
  - The terminal exon of every transcript of the gene is updated, not just one.
  - The extension report lists every gene with StringTie evidence, including rejected genes and the reason.
  - Maximum extension lowered from 10,000 bp to 2,000 bp. StringTie `-c 3 -f 0.1 -g 50`.
- **v13 (= v1.0.0): current thresholds.**
  - StringTie `-c` lowered from 3 to 2.
  - Extension minimum coverage lowered from 3.0 to 2.0.
  - Maximum extension raised from 2,000 bp to 3,000 bp.

### Current thresholds (v1.0.0)

| Step | Setting | Value |
|---|---|---|
| StringTie | min. coverage (`-c`) | 2 |
| StringTie | min. isoform fraction (`-f`) | 0.1 |
| StringTie | gap merge (`-g`) | 50 bp |
| Splice filter | min. coverage for new transcripts | 2.0 |
| Splice filter | exon-boundary tolerance | 5 bp |
| Extension | min. coverage | 2.0 |
| Extension | min. / max. extension | 10 bp / 3,000 bp |
| Extension | gap to next same-strand gene | 50 bp |
