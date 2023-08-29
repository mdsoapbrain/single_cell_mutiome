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
## Preprocessing [scripts](https://github.com/mdsoapbrain/single_cell_mutiome/blob/main/Preprocessing.R)
Single-cell multiome analysis.
Data quality control and normalization.
## hdWGCNA Pipeline [scripts](https://github.com/mdsoapbrain/single_cell_mutiome/blob/main/hdWGCNA.R)
### Creation of Metacell and Gene Networks
Construct networks using hdWGCNA.
### Module Identification
Identify modules based on network topology.
### Differential Module Eigengene Analysis
Perform differential analysis on module eigengenes.
### Functional Enrichment Analysis
Enrich modules for biological functions.
### Motif Overlap Network Analysis [scripts](https://github.com/mdsoapbrain/single_cell_mutiome/blob/main/motif_overlap_hdWGCNA.R)
Analyze the overlap of motifs within and across modules.

## Chromatin Accessibility
Analyze chromatin accessibility data.

## Downstream Analysis (Differential,Motif, Feature Linkage, LDSC)
### Differential Expression Gene Analysis [scripts](https://github.com/mdsoapbrain/single_cell_mutiome/blob/main/DEG&DAE&GeneEnrich.R)
Identify differentially expressed genes.
### Differential Accessibility Analysis [scripts](https://github.com/mdsoapbrain/single_cell_mutiome/blob/main/DEG&DAE&GeneEnrich.R)
Identify regions with differential chromatin accessibility.
## Functional and Trajectory Analysis [scripts](https://github.com/mdsoapbrain/single_cell_mutiome/blob/main/DEG&DAE&GeneEnrich.R)
- Perform pathway analysis.
- Conduct trajectory analysis for cell differentiation.
  
### Motif Enrichment
Perform motif enrichment analysis.

### Feature Linkage Analysis [scripts](https://github.com/mdsoapbrain/single_cell_mutiome/blob/main/Feature_linkage_analysis.R)
Link features from different data types (e.g., RNA and ATAC).
### LD Score Regression [scripts](https://github.com/mdsoapbrain/single_cell_mutiome/blob/main/prepare_LDSC.R)
Conduct LD Score regression for genetic correlation.
[LDSC](https://github.com/mdsoapbrain/single_cell_mutiome/blob/main/LDSC_analysis.sh)













