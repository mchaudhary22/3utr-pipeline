# UTRex method

Back to [README](../README.md) · Outputs: [outputs.md](outputs.md)

## The problem

10x Genomics 3′ single-cell RNA-seq reads come from the 3′ end of each transcript. Many plant genome annotations have short or missing 3′ UTRs, so the annotated gene often ends before the place where the reads are. Cell Ranger cannot assign reads that fall past the annotated gene end to any gene, so it discards them. These genes then look less expressed, or not expressed at all, in the count matrix.

## Steps

UTRex runs these steps for each sample (`docker/run_pipeline.sh`):

1. **GFF3 → GTF.** Convert the annotation with `gffread -T -E`, then strip `gene:` and `transcript:` prefixes from the IDs.
2. **Remove same-strand overlapping genes.** Genes that overlap on the same strand are grouped, and one gene per group is kept: a protein-coding gene if the GTF has biotype information, then the longest. All decisions go into an overlap report.
3. **Chromosome lengths.** Run `samtools faidx` on the genome FASTA.
4. **Baseline Cell Ranger run.** Run `cellranger mkref` and `cellranger count` on the original (overlap-filtered) reference. This gives the "before" matrix and the BAM file.
5. **Reference-guided StringTie.** Assemble transcripts from the BAM with the reference GTF as a guide (`-c 2 -f 0.1 -g 50`).
6. **gffcompare.** Compare the StringTie transcripts to the reference and give each one a class code.
7. **Splice-junction-aware filter.**
   - Keep transcripts tied to a reference gene (class codes `=`, `c`, `k`, `m`, `n`, `j`).
   - Drop opposite-strand, overlap and repeat codes (`o`, `s`, `x`, `r`).
   - Keep a new transcript (other codes) only if it has at least 2x coverage and at least one exon boundary within 5 bp of a known exon boundary on the same strand.
8. **Extend the 3′ ends.** For each reference gene with StringTie support, move only the end of the terminal exon of each of its transcripts. No new transcripts are added. An extension is accepted only if:
   - the StringTie transcript has at least 2x coverage;
   - the extension is between 10 bp and 3,000 bp;
   - it does not run past the chromosome end;
   - it stops 50 bp before the next gene on the same strand. If it would reach that gene, it is trimmed back. If less than 10 bp of extension is left after trimming, it is rejected.
9. **Final Cell Ranger run.** Run `cellranger mkref` and `cellranger count` on the extended reference. This gives the "after" matrix.

## Why move the gene end instead of adding transcripts?

If a new StringTie transcript is added next to the original one, the gene has two overlapping models. Cell Ranger may then treat reads that fit both as ambiguous and drop them. UTRex keeps one model per gene and only moves its 3′ end, so each read still points to one gene.

## Requirements

- Apptainer (or Singularity), on a Linux system or HPC cluster
- SLURM, if you use the provided job script
- The Docker image `mchaudhary22/3utr-pipeline-v13:latest`. The pipeline uses gffread, StringTie, gffcompare, samtools and Cell Ranger 10.0.0.

## Inputs in detail

```
/path/to/input/
├── GENOME_NAME.fa            # genome FASTA (.fa, .fasta or .fna)
├── GENOME_NAME.gff3          # annotation (.gff3 or .gff)
├── srr_list.txt              # only for list mode: one SRR ID per line
├── SRR_ID/
│   ├── SRR_ID_1.fastq.gz     # Read 1 (cell barcode + UMI)
│   ├── SRR_ID_2.fastq.gz     # Read 2 (cDNA)
│   └── SRR_ID_3.fastq.gz     # optional index read (used as I1)
└── SRR_ID_2nd_sample/
    └── ...
```

- `GENOME_NAME` is the name you pass with `-g`. UTRex also finds the files if the case differs (for example `-g tomato` finds `Tomato.gff3`).
- Uncompressed `.fastq` files also work.
- UTRex writes `GENOME_NAME.fa.fai` next to the FASTA, so the input folder must be writable.
- See `docker/srr_list.example.txt` for the list format. Blank lines and lines starting with `#` are skipped.

## Options in detail

| Flag | Meaning |
|---|---|
| `-s SRR_ID` | Run one sample |
| `-l SRR_LIST_FILE` | Run every sample in a list file (path inside the container, e.g. `/app/input/srr_list.txt`) |
| `-g GENOME_NAME` | Genome name; must match the FASTA and GFF3 file names (required) |
| `-c EXPECTED_CELLS` | Optional. Expected cell count, passed to Cell Ranger as `--expect-cells`. No default. |
| `-f FORCE_CELLS` | Optional. Exact cell count, passed to Cell Ranger as `--force-cells`. Overrides `-c`. |
| `-h` | Show help |

You must give either `-s` or `-l`. If you give neither `-c` nor `-f`, Cell Ranger detects the cell count itself.

## List mode

```bash
apptainer run \
  --bind /path/to/input:/app/input \
  --bind /path/to/output:/app/output \
  3utr_pipeline_v13.sif \
  /app/run_pipeline.sh -l /app/input/srr_list.txt -g GENOME_NAME -c 10000
```

Samples run one after another. The reference GTF and the baseline Cell Ranger reference are built once and reused.

## SLURM

Edit `script/run_utrex_slurm.sh` (allocation, email, paths, `SRR_ID`, `GENOME_NAME`), then:

```bash
sbatch script/run_utrex_slurm.sh
```

The script requests 16 CPUs, 128 GB of memory and 48 hours; change these for your data and cluster.

## Re-running

Each step is skipped if its output already exists. To redo a step, delete its output first.
