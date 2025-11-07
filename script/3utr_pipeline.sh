#!/bin/bash
#SBATCH --job-name=v1_rice
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --output=v1_rice_SRR13853439_%j.out
#SBATCH --error=v1_rice_SRR13853439_%j.err
#SBATCH --mail-user=mchaudhary@vt.edu
#SBATCH --mail-type=ALL

# Load Apptainer module
module load apptainer/1.4.0

# Set working directory
cd /projects/songli_lab/Manisha_3UTR_2025/3utr_docker

# Define input/output paths
INPUT_DIR="/projects/songli_lab/Manisha_3UTR_2025/3utr_docker/testing/input"
OUTPUT_DIR="/projects/songli_lab/Manisha_3UTR_2025/3utr_docker/output/rice_v1"

# Check if .sif image already exists, otherwise build from new Docker image
if [[ ! -f "3utr_pipeline_v1.sif" ]]; then
    echo "Building new Apptainer image from Docker Hub..."
    apptainer build 3utr_pipeline_v1.sif docker://mchaudhary22/3utr-pipeline-v1:latest
fi

# Run the pipeline for a single sample
apptainer run \
  --bind ${INPUT_DIR}:/app/input \
  --bind ${OUTPUT_DIR}:/app/output \
  3utr_pipeline_v1.sif \
  /app/run_pipeline.sh -s SRR13853439 -c 1800000 -g rice
