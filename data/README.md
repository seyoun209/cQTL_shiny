# Data Directory

This directory should contain the following data files that are not tracked in Git due to their large size:

- `counts_annot.Rdata` 
- `expression_data.json`
- `genes.json`
- `normalized_log2_gene.txt`
- `raw_counts_with_coords.csv`

You can generate these files by running the script in `scripts/00_cleanup_data.R`.
