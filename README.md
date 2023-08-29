# single_cell_mutiome: Single-Cell RNA-Seq and ATAC-Seq Analysis Pipeline

## Overview
This pipeline is designed for the integrated analysis of single-cell RNA-Seq (scRNA) and single-cell ATAC-Seq (scATAC) data. It covers everything from data demultiplexing to downstream analyses, including network construction, differential expression, and functional enrichment.

![My Image](https://github.com/mdsoapbrain/single_cell_mutiome/blob/main/scmutiome_pipeline.png)


## Table of Contents
- Data Preparation
- Preprocessing
- hdWGCNA Pipeline
- Chromatin Accessibility
- Gene Expression
- Differential Analyses
- Functional and Trajectory Analysis
- Motif Enrichment


## Data Preparation
scATAC (fastq.gz) and scRNA (fastq.gz)
Demultiplex raw sequencing data using Cell Ranger - Arc.
Cell Ranger Phase
Make a reference with mRatBN7.
Obtain counts using Cell Ranger count.
## Preprocessing
Single-cell multiome analysis.
Data quality control and normalization.
## hdWGCNA Pipeline
### Creation of Metacell and Gene Networks
Construct networks using hdWGCNA.
### Module Identification
Identify modules based on network topology.
### Differential Module Eigengene Analysis
Perform differential analysis on module eigengenes.
### Functional Enrichment Analysis
Enrich modules for biological functions.
### Motif Overlap Network Analysis
Analyze the overlap of motifs within and across modules.

## Chromatin Accessibility
Analyze chromatin accessibility data.


## Differential Analyses
### Differential Expression Gene Analysis
Identify differentially expressed genes.
### Differential Accessibility Analysis
Identify regions with differential chromatin accessibility.
## Functional and Trajectory Analysis
- Perform pathway analysis.
- Conduct trajectory analysis for cell differentiation.
  
## Motif Enrichment
Perform motif enrichment analysis.

## Downstream Analysis
### Feature Linkage Analysis
Link features from different data types (e.g., RNA and ATAC).
### LD Score Regression
Conduct LD Score regression for genetic correlation.














