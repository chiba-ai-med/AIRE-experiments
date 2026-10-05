# Dockerfile for AIRE-experiments
# Multi-stage build, mirrors /home/koki/dev/eCCI-experiments style.
# Used with `snakemake --use-singularity` once the image is published as
# koki/aire_workflow:<tag> (matches DOCKER_IMAGE in workflow/preprocess.smk).

FROM condaforge/mambaforge:latest AS builder

# --- R env (mirrors workflow/envs/r-machima.yaml) ---
RUN mamba install -y -c conda-forge -c bioconda \
    r-base=4.4 \
    r-matrix r-dplyr r-tidyr r-readr \
    r-ggplot2 r-viridis \
    r-remotes r-biocmanager \
    && mamba clean -afy

# Heavy R packages (CRAN/Bioc/GitHub) installed via setup script.
COPY setup_r_packages.R /tmp/setup_r_packages.R
RUN Rscript /tmp/setup_r_packages.R

# --- Python env (mirrors workflow/envs/py-hic.yaml) ---
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.11 \
    numpy pandas scipy h5py pyranges \
    scanpy anndata mudata macs2 \
    cooler cooltools pairtools \
    wget aria2 \
    && mamba clean -afy
RUN pip install --no-cache-dir snapatac2

#################################
FROM condaforge/mambaforge:latest

COPY --from=builder /opt/conda /opt/conda
WORKDIR /work
CMD ["/bin/bash"]

LABEL maintainer="AIRE-experiments"
LABEL description="Cell-type deconvolution validation (Machima2) on multiome + Hi-C"
LABEL version="0.1.0"
