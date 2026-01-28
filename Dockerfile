FROM continuumio/miniconda3:latest

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
# We include snakemake inside so the container is self-orchestrating if needed
RUN conda install -y \
    snakemake \
    bcftools \
    pysam \
    tabix \
    samtools \
    pypgx \
    pandas \
    openjdk=11  # Required for Picard

# Create a working directory
WORKDIR /pipeline

# Copy pipeline files (optional, depends on how you run it)
# COPY . /pipeline

CMD ["/bin/bash"]
