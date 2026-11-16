#!/bin/bash
#SBATCH --job-name=utrex_<SRR_ID>
#SBATCH --account=<YOUR_ALLOCATION>
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --output=utrex_<SRR_ID>_%j.out
#SBATCH --error=utrex_<SRR_ID>_%j.err
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL

# Edit everything in <...>, every /path/to/..., SRR_ID and GENOME_NAME before submitting.
# Partition name and module version depend on your cluster.

module load apptainer/1.4.0
cd /path/to/workdir



INPUT_DIR="/path/to/input"
OUTPUT_DIR="/path/to/output"
SRR_ID="SRR0000001"          # folder name under INPUT_DIR holding the FASTQs
GENOME_NAME="your_genome"    # matches your_genome.fa and your_genome.gff3 in INPUT_DIR

mkdir -p "$OUTPUT_DIR"



if [[ ! -f "3utr_pipeline_v13.sif" ]]; then
    apptainer pull 3utr_pipeline_v13.sif docker://mchaudhary22/3utr-pipeline-v13:latest
fi


apptainer run \
  --bind ${INPUT_DIR}:/app/input \
  --bind ${OUTPUT_DIR}:/app/output \
  3utr_pipeline_v13.sif \
  /app/run_pipeline.sh -s "$SRR_ID" -g "$GENOME_NAME"
