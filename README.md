# Chondrocyte Gene Expression Explorer

This repository contains tools for exploring gene expression data in chondrocytes.

## Quick Start

1. **Static HTML version**:
   - Open `www/interactive.html` in your browser

2. **Shiny app version**:
   - Run `./run_app.sh`
   - Open http://localhost:8888 in your browser

## Data Files

- `data/counts_annot.Rdata`: Annotated expression data
- `data/genes.json`: List of significant genes
- `data/expression_data.json`: Expression data for visualization

## Scripts

- `scripts/00_cleanup_data.R`: Data processing script
- `app.R`: Shiny application

