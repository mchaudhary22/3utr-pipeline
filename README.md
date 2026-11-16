# UTRex: 3′ UTR extension for plant single-cell RNA-seq references

UTRex uses your 10x single-cell reads to extend short or missing 3′ UTRs in a plant genome annotation, so Cell Ranger counts reads it would otherwise discard.

## How it works

- Clean the annotation: remove genes that overlap on the same strand.
- Run Cell Ranger, then StringTie, to find where the reads actually end.
- Filter out noisy StringTie transcripts.
- Move only the 3′ end of each supported gene (no new transcripts are added).
- Re-run Cell Ranger on the extended reference.

## Quick start

```bash
apptainer pull 3utr_pipeline_v13.sif docker://mchaudhary22/3utr-pipeline-v13:latest

apptainer run \
  --bind /path/to/input:/app/input \
  --bind /path/to/output:/app/output \
  3utr_pipeline_v13.sif \
  /app/run_pipeline.sh -s SRR0000001 -g GENOME_NAME
```

For HPC, edit and submit `script/run_utrex_slurm.sh`.

## Inputs

```
/path/to/input/
├── GENOME_NAME.fa
├── GENOME_NAME.gff3
└── SRR_ID/
    ├── SRR_ID_1.fastq.gz
    ├── SRR_ID_2.fastq.gz
    └── SRR_ID_3.fastq.gz   # optional
```

## Main outputs

- `SRR_ID_matrix/`: Cell Ranger result before extension
- `SRR_ID_stringtie_matrix/`: Cell Ranger result after extension
- `SRR_ID_extension_report.tsv`: which genes were extended, by how much, or why not

## Options

| Flag | Meaning |
|---|---|
| `-s SRR_ID` | One sample |
| `-l LIST_FILE` | Several samples, one SRR ID per line |
| `-g GENOME_NAME` | Genome name, matching the FASTA and GFF3 (required) |
| `-c N` | Expected cells (optional, no default) |
| `-f N` | Force exact cell count (optional) |
| `-h` | Help |

## Cell Ranger

UTRex uses Cell Ranger 10.0.0 (10x Genomics software under its own license). Details on obtaining it will be added.

## More

- Full method: [docs/method.md](docs/method.md)
- All outputs and report columns: [docs/outputs.md](docs/outputs.md)
- Changelog: [CHANGELOG.md](CHANGELOG.md)

## Citation

Manuscript in preparation.
