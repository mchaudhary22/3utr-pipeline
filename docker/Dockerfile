# Use Ubuntu as base image
FROM ubuntu:20.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/conda/bin:$PATH"

# Install system dependencies (added gawk)
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    python3 \
    python3-pip \
    unzip \
    gzip \
    tar \
    dos2unix \
    gawk \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

# Initialize conda and accept TOS
RUN conda init bash && \
    conda config --set channel_priority strict && \
    echo "yes" | conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    echo "yes" | conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install bioinformatics tools via conda
RUN conda install -c conda-forge -c bioconda \
    gffread \
    stringtie \
    gffcompare \
    sra-tools \
    samtools \
    && conda clean -a

# Install Cell Ranger
# Cell Ranger is not redistributed with UTRex (10x Genomics license).
# Download Cell Ranger 10.0.0 yourself from:
#   https://www.10xgenomics.com/support/software/cell-ranger/downloads
# and place the tarball next to this Dockerfile with this exact filename:
#   cellranger-10.0.0.tar.gz
WORKDIR /opt
COPY cellranger-10.0.0.tar.gz /opt/
RUN tar -xzf cellranger-10.0.0.tar.gz && \
    rm cellranger-10.0.0.tar.gz

# Add Cell Ranger to PATH
ENV PATH="/opt/cellranger-10.0.0:$PATH"

# Create working directory
WORKDIR /app

# Copy pipeline scripts
COPY run_pipeline.sh /app/
COPY extend_utr_boundaries.py /app/
COPY remove_overlap.py /app/
COPY filter_stringtie_splicing.py /app/
COPY requirements.txt /app/

# Convert Windows line endings to Unix and make scripts executable
RUN dos2unix /app/*.sh /app/*.py && \
    chmod +x /app/*.sh /app/*.py

# Install Python requirements
RUN pip3 install -r requirements.txt

# Create directories for input/output
RUN mkdir -p /app/input /app/output /app/reference

# Set default command
CMD ["/app/run_pipeline.sh"]