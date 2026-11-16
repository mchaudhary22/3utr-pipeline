# UTRex outputs

Back to [README](../README.md) · Method: [method.md](method.md)

All paths are inside the output folder (`/app/output`). `GENOME_NAME` and `SRR_ID` are your values.

## Shared reference files (`reference/`)

- `reference/GENOME_NAME.gtf`: reference GTF after prefix stripping and overlap removal
- `reference/GENOME_NAME_overlap_report.tsv`: overlap report: every gene, with `kept` or `removed` and which gene it was kept in favor of
- `reference/GENOME_NAME_genome/`: baseline Cell Ranger reference
- `reference/SRR_ID_stringtie_genome/`: extended Cell Ranger reference for this sample
- `reference/mkref_*/`: Cell Ranger mkref working folders (created by Cell Ranger)

## Per-sample files

- `SRR_ID/`: uncompressed FASTQs renamed for Cell Ranger (`SRR_ID_S1_L001_R1_001.fastq`, `R2`, `I1`)
- `SRR_ID_matrix/`: baseline Cell Ranger count ("before" matrix and `outs/possorted_genome_bam.bam`)
- `SRR_ID_stringtie.gtf`: StringTie assembly
- `SRR_ID_stringtie_srr_annotation.*`: gffcompare output (`.annotated.gtf`, `.loci`, `.stats`, `.tracking`, `.refmap`, `.tmap`)
- `SRR_ID_stringtie_filtered.gtf`: StringTie transcripts that passed the splice filter
- `SRR_ID_splice_filter_report.tsv`: splice filter report: one row per transcript with `class_code`, `coverage`, `has_splice`, `decision` (`kept`/`removed`) and `reason`
- `SRR_ID_reference_extended.gtf`: reference GTF with extended 3′ ends
- `SRR_ID_extension_report.tsv`: extension report (see below)
- `SRR_ID_stringtie_matrix/`: final Cell Ranger count on the extended reference ("after" matrix)

## In the input folder

- `GENOME_NAME.fa.fai`: FASTA index from `samtools faidx`

## Extension report

`SRR_ID_extension_report.tsv` has one row per reference gene that had StringTie support.

| Column | Meaning |
|---|---|
| `gene_id` | Reference gene ID |
| `strand` | `+` or `-` |
| `old_3prime` | Annotated 3′ coordinate (gene end on `+`, gene start on `-`) |
| `new_3prime` | 3′ coordinate after extension |
| `extension_bp` | Length of the extension in bp |
| `coverage` | StringTie coverage of the transcript used |
| `status` | `extended`, `extended_trimmed` or `rejected` |
| `reason` | Why (see below) |

- **`extended`**: the gene end was moved to the StringTie end. `reason` is `ok`.
- **`extended_trimmed`**: the StringTie end would have come within 50 bp of the next gene on the same strand, so the extension was shortened to stop 50 bp before that gene. `reason` is `neighbor_cap_applied`.
- **`rejected`**: the gene was not changed. `old_3prime`, `new_3prime` and `extension_bp` are `.`. `reason` is one of:
  - `too_long_<N>bp`: extension over 3,000 bp
  - `too_short_<10bp`: extension under 10 bp
  - `wrong_direction`: StringTie does not reach past the annotated 3′ end
  - `exceeds_chrom_length_<N>` / `exceeds_chrom_start`: would run off the chromosome
  - `neighbor_overlap_cap_<N>`: less than 10 bp left after trimming for the neighbor gene
  - `low_coverage_<N>x`: coverage under 2x
  - `no_ref_boundary`: gene not found in the reference GTF
