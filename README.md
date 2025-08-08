# 3utr-pipeline
A one-step 3' UTR single-cell RNA-seq analysis pipeline using Cell Ranger and StringTie.

## Overview
This pipeline performs the following steps:

1. **Initial Cell Ranger count** - Standard gene expression quantification
2. **CDS-only analysis** - Extracts coding sequence regions and creates CDS-specific reference
3. **StringTie assembly** - Identifies novel transcript isoforms and 3' UTR regions
4. **Final quantification** - Re-quantifies using the improved reference with identified 3' UTRs

## Docker Image
Pre-built Docker image is available at: 
```
mchaudhary22/3utr-pipeline:latest
```

## Input Requirements
Your input directory must contain:

### For each sample:
**Sequencing data**: Place SRR data in folders named with SRR IDs (e.g., SRR8257100/)
    
- SRR_ID_1.fastq.gz (or .fastq) - Read 1
- SRR_ID_2.fastq.gz (or .fastq) - Read 2
- SRR_ID_3.fastq.gz (or .fastq) - Index reads (optional)

### Reference files:
- **Genome FASTA**: GENOME_NAME.fa (or .fasta)
- **Annotation GFF3**: GENOME_NAME.gff3 (or .gff)

###  For multiple samples:
**Sample list**: srr_list.txt containing one SRR ID per line

## Output

The pipeline generates several outputs for each sample:

- {SRR_ID}_matrix/ - Initial Cell Ranger results
- {SRR_ID}_cds2exon_matrix/ - CDS-only analysis results
- {SRR_ID}_stringtie.gtf - StringTie assembled transcripts
- {SRR_ID}_stringtie_matrix/ - Final results with 3' UTR quantification

## Usage
  
### Prerequisites
- Apptainer/Singularity installed on your system
- Input files prepared as described above

### Running with Apptainer
#### Single Sample
```bash

# Build the container
apptainer build 3utr_pipeline.sif docker://mchaudhary22/3utr-pipeline:latest
      
# Run pipeline
apptainer run \
  --bind /path/to/input:/app/input \
  --bind /path/to/output:/app/output \
  3utr_pipeline.sif \
  /app/run_pipeline.sh -s <SRR_Sample_ID> -c <Expected_Cells> -g <Genome_Name>

     
#### Multiple Samples
```bash
apptainer run \
  --bind /path/to/input:/app/input \
  --bind /path/to/output:/app/output \
  3utr_pipeline.sif \
  /app/run_pipeline.sh -l /app/input/srr_list.txt -c <Expected_Cells> -g <Genome_Name> 


### Parameters
-s SRR_ID: Single SRR accession to process
-l SRR_LIST_FILE: Path to file containing multiple SRR IDs (one per line)
-c EXPECTED_CELLS: Expected number of cells (default: 1800000)
-g GENOME_NAME: Species/genome name (must match your input file names)
-h: Display help message


## SLURM Integration
For HPC environments, use the provided SLURM script and modify the paths and parameters according to your system:
```bash
# Edit the script paths and parameters
sbatch script/3utr_pipeline.sh


