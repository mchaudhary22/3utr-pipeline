#!/bin/bash
#SBATCH --job-name=3utr_pipeline
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --output=3utr_pipeline_%j.out
#SBATCH --error=3utr_pipeline_%j.err
#SBATCH --mail-user=mchaudhary@vt.edu
#SBATCH --mail-type=ALL

# Load apptainer module
module load apptainer/1.4.0

# Set working directory
cd /projects/songli_lab/Manisha_Canker_2025/3UTR_Project/3utr_docker/testing

# Define paths
INPUT_DIR="/projects/songli_lab/Manisha_Canker_2025/3UTR_Project/3utr_docker/testing/input"
OUTPUT_DIR="/projects/songli_lab/Manisha_Canker_2025/3UTR_Project/3utr_docker/testing/output"

# Pull and convert Docker image to Apptainer
if [[ ! -f "3utr_pipeline.sif" ]]; then
    apptainer build 3utr_pipeline.sif docker://mchaudhary22/3utr-pipeline:latest
fi

# Run the pipeline
apptainer run \
  --bind ${INPUT_DIR}:/app/input \
  --bind ${OUTPUT_DIR}:/app/output \
  3utr_pipeline.sif \
  /app/run_pipeline.sh -s SRR13853439 -c 1800000 -g rice

# For multiple samples, use:
# apptainer run \
#   --bind ${INPUT_DIR}:/app/input \
#   --bind ${OUTPUT_DIR}:/app/output \
#   3utr_pipeline.sif \
#   /app/run_pipeline.sh -l /app/input/srr_list.txt -c 1800000 -g rice
