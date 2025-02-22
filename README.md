# T Cell Activation ATAC-seq Visualization App
![Docker](https://img.shields.io/badge/docker-ready-blue.svg)
![R-Shiny](https://img.shields.io/badge/R--Shiny-app-brightgreen.svg)
![License](https://img.shields.io/badge/license-MIT-blue.svg)

An interactive web application for visualizing and analyzing ATAC-seq data from T cell activation studies. This application provides insights into chromatin accessibility changes during T cell activation and immune responses.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Data Description](#data-description)
- [Docker Deployment](#docker-deployment)
- [Usage Guide](#usage-guide)
- [Configuration](#configuration)
- [Contributing](#contributing)

## Overview
This Shiny application enables immunologists and researchers to explore chromatin accessibility changes during T cell activation. It integrates with processed data from the [ATACseq_Pipeline](../ATACseq_Pipeline) to provide interactive visualization and analysis capabilities focused on immune cell dynamics.

## Features
- Interactive visualization of T cell ATAC-seq data
- Time-course analysis of chromatin accessibility
- Integration with T cell activation markers
- Dynamic comparison of activation states
- Statistical analysis tools
- Custom region analysis
- Data export functionality
- Publication-ready visualizations

## Quick Start

Using Docker:
```bash
docker pull rohitrrj/tcell-atac-viz
docker run -p 3838:3838 rohitrrj/tcell-atac-viz
```

Access the application at: http://localhost:3838

## Installation

### Local Installation
1. Clone the repository:
```bash
git clone https://github.com/rohitrrj/t_cell_activation_atac.git
cd t_cell_activation_atac
```

2. Install R dependencies:
```R
install.packages(c(
  "shiny",
  "ggplot2",
  "DT",
  "plotly",
  "dplyr",
  "tidyr",
  "GenomicRanges"
))
```

3. Run the application:
```R
shiny::runApp("app")
```

## Data Description

### Input Data Types
- Processed ATAC-seq peak files
- Time-course accessibility data
- T cell activation markers
- Gene expression correlation data
- Pathway annotations

### Data Structure
```
app/
├── data/
│   ├── peak_data/          # ATAC-seq peaks
│   ├── activation_markers/ # T cell markers
│   ├── time_course/       # Time series data
│   └── annotations/       # Genomic annotations
└── www/
    └── documentation/
```

## Docker Deployment

### Building the Image
```dockerfile
FROM rocker/shiny:4.1.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev

# Install Bioconductor packages
RUN R -e 'install.packages("BiocManager"); \
    BiocManager::install(c("GenomicRanges", "rtracklayer"))'

# Copy application files
COPY app /srv/shiny-server/app
COPY shiny-server.sh /usr/bin/shiny-server.sh
COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

# Set permissions
RUN chmod +x /usr/bin/shiny-server.sh

# Expose port
EXPOSE 3838

# Start Shiny server
CMD ["/usr/bin/shiny-server.sh"]
```

### Running the Container
```bash
# Build image
docker build -t tcell-atac-viz .

# Run container
docker run -d \
  -p 3838:3838 \
  --name tcell-atac-viz \
  tcell-atac-viz
```

## Usage Guide

1. **Data Exploration**
   - Select time points
   - Choose cell populations
   - Compare activation states
   - Filter regions of interest

2. **Visualization Options**
   - Time-course plots
   - Peak accessibility heatmaps
   - Genomic region browser
   - Pathway enrichment visualization

3. **Analysis Features**
   - Differential accessibility testing
   - Time-course clustering
   - Motif enrichment analysis
   - Integration with gene expression

4. **Export Options**
   - High-resolution plots
   - Data tables
   - Analysis reports
   - Region coordinates

## Configuration

### Application Settings
Modify `app/app.R` for:
```R
# App configuration
options(
  shiny.maxRequestSize = 200*1024^2,
  shiny.plot.width = 1000,
  shiny.plot.height = 800
)

# Analysis parameters
SIGNIFICANCE_THRESHOLD = 0.05
MIN_FOLD_CHANGE = 2.0
TIME_POINTS = c(0, 2, 6, 24, 48)
```

### Server Configuration
Edit `shiny-customized.config`:
```
run_as shiny;
preserve_logs true;
access_log /var/log/shiny-server/access.log tiny;
server {
  listen 3838;
  location / {
    site_dir /srv/shiny-server/app;
    log_dir /var/log/shiny-server;
    directory_index on;
  }
}
```

## Analysis Workflows

### Time-Course Analysis
1. Select activation time points
2. Choose normalization method
3. Apply clustering algorithm
4. Visualize temporal patterns

### Differential Accessibility
1. Define comparison groups
2. Set significance thresholds
3. Identify differential regions
4. Perform pathway analysis

### Region Analysis
1. Upload regions of interest
2. Compare across conditions
3. Calculate enrichment scores
4. Generate summary statistics

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

### Development Guidelines
- Follow R code style guidelines
- Add unit tests for new features
- Update documentation
- Test Docker deployment

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Applications
This visualization tool has been used in the following publications:

1. "PD-1 combination therapy with IL-2 modifies CD8+ T cell exhaustion program"
   - *Nature*. 2022 Oct;610(7933):737-743
   - DOI: [10.1038/s41586-022-05257-0](https://doi.org/10.1038/s41586-022-05257-0)
   - PMID: [36215562](https://pubmed.ncbi.nlm.nih.gov/36215562)
   - PMCID: [PMC9927214](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9927214/)
   - Used for analyzing T cell activation dynamics and chromatin states

2. "Epigenetic signature of PD-1+ TCF1+ CD8 T cells that act as resource cells during chronic viral infection"
   - *Proc Natl Acad Sci U S A*. 2022 Feb 22;119(8):e2117314119
   - DOI: [10.1073/pnas.2117314119](https://doi.org/10.1073/pnas.2117314119)
   - PMID: [35085847](https://pubmed.ncbi.nlm.nih.gov/35085847)
   - Applied for studying chromatin accessibility changes during T cell responses

Code availability:
⭐ [rohitrrj/t_cell_activation_atac](https://github.com/rohitrrj/t_cell_activation_atac) - Interactive visualization tool for T cell activation ATAC-seq analysis

## Related Projects
- [ATACseq_Pipeline](../ATACseq_Pipeline)
- [macrophage_atac](../macrophage_atac)
- [rheumatoid_arthritis_atac](../rheumatoid_arthritis_atac)

## Acknowledgments
- Built with R Shiny framework
- Data processing: ATACseq_Pipeline
- Supporting institutions and funding
- T cell biology expertise and guidance
