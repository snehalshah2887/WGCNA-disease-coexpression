# =============================================================================
# Dockerfile — 22q13 WGCNA Pipeline
#
# Base: Bioconductor 3.18 / R 4.3
# WGCNA is no longer available on CRAN (deprecated); it is installed via
# Bioconductor where it is actively maintained.
#
# Build:
#   docker build -t 22q13-wgcna:latest .
#
# Run a single module (example):
#   docker run --rm \
#     -v $(pwd)/data/input:/input:ro \
#     -v $(pwd)/results:/results \
#     22q13-wgcna:latest \
#     Rscript /pipeline/scripts/01_preprocess.R \
#       /input/gene_expression_filtered.csv \
#       /results/01_preprocessed
# =============================================================================

FROM bioconductor/bioconductor_docker:RELEASE_3_18

LABEL maintainer="Snehal Shah <ecosnehal@gmail.com>"
LABEL description="22q13 / Phelan-McDermid Syndrome WGCNA analysis pipeline"
LABEL org.opencontainers.image.source="https://github.com/YOUR_USERNAME/22q13_WGCNA"

# System dependencies for R packages and headless graphics
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    xvfb \
    && rm -rf /var/lib/apt/lists/*

# Set CRAN mirror
ENV R_CRAN_MIRROR=https://cloud.r-project.org

# Install CRAN packages
RUN Rscript -e "install.packages(c( \
    'dplyr', \
    'openxlsx', \
    'readtext', \
    'tidyr', \
    'ggplot2', \
    'gplots', \
    'RColorBrewer', \
    'VennDiagram', \
    'venn', \
    'stringr', \
    'forcats', \
    'matrixStats', \
    'remotes' \
  ), repos = Sys.getenv('R_CRAN_MIRROR'), ask = FALSE)"

# Install Bioconductor packages
RUN Rscript -e "BiocManager::install(c( \
    'WGCNA', \
    'RCy3', \
    'biomaRt', \
    'biomartr', \
    'org.Hs.eg.db', \
    'preprocessCore', \
    'impute', \
    'GO.db', \
    'AnnotationDbi' \
  ), ask = FALSE, update = FALSE)"

# Copy pipeline scripts into the image
COPY scripts/ /pipeline/scripts/

# Use Xvfb for headless display (needed for any interactive graphic calls)
ENV DISPLAY=:99
RUN echo '#!/bin/bash\nXvfb :99 -screen 0 1024x768x24 &\nexec "$@"' > /entrypoint.sh \
    && chmod +x /entrypoint.sh

WORKDIR /pipeline
ENTRYPOINT ["/entrypoint.sh"]
CMD ["bash"]
