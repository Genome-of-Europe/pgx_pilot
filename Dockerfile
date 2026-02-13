FROM continuumio/miniconda3:latest
LABEL org.opencontainers.image.source=https://github.com/Genome-of-Europe/pgx_pilot

# Set non-interactive mode
ENV DEBIAN_FRONTEND=noninteractive

# Install system basics
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Setup Conda Channels
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Install Bioinformatics Tools & Snakemake
RUN conda install -y \
    python=3.11 \
    snakemake \
    bcftools \
    pysam \
    tabix \
    samtools \
    pypgx \
    pandas \
    openjdk=11

# Create a working directory
WORKDIR /pipeline

# Copy pipeline files into the image
COPY . /pipeline

# Ensure all files are readable and executable (for Singularity)
RUN chmod -R a+rX /pipeline

# Default command
CMD ["/bin/bash"]
