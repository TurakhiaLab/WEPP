# Start from the official Ubuntu 24.04.2 LTS base image
FROM ubuntu:24.04

# Set environment variables to avoid interactive prompts during install
ENV DEBIAN_FRONTEND=noninteractive

# Install essential packages
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    vim \
    sudo \
    pip \
    build-essential \
    libopenmpi-dev \
    libboost-all-dev \
    python3-pandas \
    pkg-config \
    zip \
    cmake \
    libprotobuf-dev \
    protobuf-compiler \
    libtbb-dev \
    python3-pip \
    git \
    nodejs \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Download and install Miniforge (Conda)
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" \
    && bash Miniforge3-$(uname)-$(uname -m).sh -b \
    && rm Miniforge3-$(uname)-$(uname -m).sh

# Set Conda/Miniforge paths
ENV PATH=/root/miniforge3/bin:$PATH

# Initialize Conda and ensure it’s sourced on login
RUN /root/miniforge3/bin/conda init bash && \
    echo "source ~/.bashrc" >> /root/.bash_profile

# Install Snakemake using pip
RUN pip install snakemake

RUN echo 'run-wepp() {' >> /root/.bashrc && \
    echo '    /root/miniforge3/bin/snakemake -s /workspace/workflow/Snakefile "$@"' >> /root/.bashrc && \
    echo '}' >> /root/.bashrc && \
    echo 'export -f run-wepp' >> /root/.bashrc

# Default command: start an interactive login shell so .bashrc is sourced
CMD [ "bash", "--login" ]
