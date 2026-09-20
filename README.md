# scDiagnostics Interactive App

A Shiny web application for single-cell RNA sequencing annotation diagnostics.

## Citation

If you use the `scDiagnostics` app in published research, please cite:

Christidis A, Ghazi A, Chawla S, Turaga N, Gentleman R, Geistlinger L (2026). scDiagnostics: systematic assessment of cell type annotation in single-cell transcriptomics data. *Briefings in Bioinformatics*, 27(5), bbag496. doi: [10.1093/bib/bbag496](https://doi.org/10.1093/bib/bbag496).

## App Link

You may access the Shiny application here: [scDiagnostics App Link](https://ccb.connect.hms.harvard.edu/scDiagnosticsApp/)

## What is this?

This app provides an interactive interface for the `scDiagnostics` R/Bioconductor package, allowing users to assess the quality of single-cell annotation through a web browser.

## Features

- **PCA Projections**: Visualize cell type distributions in principal component space
- **Discriminant Projections**: Visualize cell type distributions in discriminant space
- **Anomaly Detection**: Detect annotation inconsistencies using isolation forests
- **Graph Integration**: Detect annotation inconsistencies using graph analysis
- **Interactive Plots**: Real-time parameter adjustment and visualization
- **Download Results**: Save high-resolution plots

## Usage

1. Select your datasets in the Data Overview tab
2. Configure parameters and diagnostic plots and statistical measures
3. Download results as needed