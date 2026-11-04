# Use Miniforge as modern, maintained base with mamba pre-configured
FROM condaforge/miniforge3:latest

LABEL org.opencontainers.image.source=https://github.com/Genome-of-Europe/pgx_pilot

# Metadata and non-interactive settings
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install essential system libraries in one layer and clean apt cache
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    git \
    wget \
    curl \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set working directory early to organize the build context
WORKDIR /pipeline

# Pre-create directory scaffolding for runtime mounts
RUN mkdir -p /pipeline/resources /pipeline/data /pipeline/results/temp

# Copy ONLY the environment file first to leverage Docker layer caching
COPY env.yml /tmp/env.yml

# Create the environment, clean up mamba cache
RUN mamba env create -f /tmp/env.yml -n pgx_pilot && mamba clean -afy

# Ensure the environment is on the PATH
ENV PATH=/opt/conda/envs/pgx_pilot/bin:$PATH

# Copy the rest of the pipeline code
COPY . /pipeline

# Permissions for HPC/Singularity compatibility
RUN chmod -R a+rX /pipeline

CMD ["/bin/bash"]
