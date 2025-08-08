# Use Ubuntu as base image
FROM ubuntu:20.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/conda/bin:$PATH"

# Install system dependencies
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
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

# Initialize conda
RUN conda init bash

# Install bioinformatics tools via conda
RUN conda install -c conda-forge -c bioconda \
    gffread \
    stringtie \
    sra-tools \
    && conda clean -a

# Download and install Cell Ranger
WORKDIR /opt
RUN wget -O cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1752564662&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=U9r7pOLdfiTHkQotlUBAn2XAh6lga162Bekme5lTbzH6HoTrJnVh4eJOq8zJ9fID8Ilzx2LRjNhPw-OXuMYX2hAtmsgQBvETgBoC2ngek~-~Txi-HPW7nvxRfsbPpbpNuwSpk31jFz8rmszcHccb~TlACNkQLpdN-mKDLoRYj9xklMQwhbSYhXqkXt2wYFIV67rXIIAyfxJp-qx0gVKErhej~gY7G9qpgY6Naz5X7hx~HoGU1qZuFA1HLPlkbO58Ck2LV8WM2KPCCo0WDxfDnizbs1KsUOe4Gl2EpI8gJMFxelXuaeVgFOBTjqoZpolVDAjAVAQGeYn0~r66QgAPAg__" && \
    tar -xzf cellranger-9.0.1.tar.gz && \
    rm cellranger-9.0.1.tar.gz

# Add Cell Ranger to PATH
ENV PATH="/opt/cellranger-9.0.1:$PATH"

# Create working directory
WORKDIR /app

# Copy pipeline scripts
COPY run_pipeline.sh /app/
COPY update_ref_gene_format.py /app/
COPY requirements.txt /app/

# Convert Windows line endings to Unix and make scripts executable
RUN dos2unix /app/*.sh /app/*.py && \
    chmod +x /app/*.sh

# Install Python requirements
RUN pip3 install -r requirements.txt

# Create directories for input/output
RUN mkdir -p /app/input /app/output /app/reference

# Set default command
CMD ["/app/run_pipeline.sh"]