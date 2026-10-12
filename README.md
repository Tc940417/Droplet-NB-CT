# Droplet NB-CT

**Droplet NB-CT: A streamlined and highly versatile single-cell chromatin immunoprecipitation method for plants**

## Overview

Droplet NB-CT is a single-cell histone modification profiling method developed for plant tissues. By incorporating nanoscale concanavalin A (ConA) magnetic beads into the Nano-CT framework, this approach enables efficient nuclear capture and high-quality single-cell chromatin profiling.

This repository contains the computational scripts used for processing and analyzing **DNBelab C4 (DNBC4) sequencing data**, with a focus on two histone modifications, **H3K4me3 and H3K27ac**, in rice (*Oryza sativa*) and Arabidopsis (*Arabidopsis thaliana*).

The analyses include sequencing quality assessment, peak calling, cell-type annotation, chromatin profiling, motif activity analysis, multi-omics integration, and identification of cell-type-specific epigenomic responses to heat stress.

## Repository Structure

```text
Droplet-NB-CT/
├── DNBC4/                   # DNBC4 sequencing preprocessing
├── Assay_pre/               # Peak-by-cell matrices and chromVAR
├── MACS/                    # Histone modification peak calling
├── Bulk_chip/               # Bulk ChIP-seq preprocessing
├── DeepTools/               # Signal analysis and visualization
├── Sequencing saturation/   # Sequencing saturation analysis
├── sctype/                  # Cell-type annotation
├── GLUE/                    # Multi-omics integration
├── Cicero/                  # Chromatin regulatory analysis
├── fig1/                    # Quality control and benchmarking
├── fig2/                    # Cell-type identification
├── fig3/                    # Cell-type-specific chromatin features
├── fig4/                    # Heat-stress response analysis
├── README.md
└── LICENSE
```
# Software Versions

This document lists the software and package versions used for the computational analyses associated with the **Droplet NB-CT** project.

## 1. R Environment

| Software | Version |
|---|---|
| R | 4.4.2 (2024-10-31 ucrt) |

### Single-cell and Epigenomic Analysis

| R Package | Version |
|---|---|
| Seurat | 5.5.0 |
| SeuratObject | 5.4.0 |
| Signac | 1.17.1 |
| SCP | 0.5.6 |
| MAST | 1.32.0 |
| DESeq2 | 1.46.0 |
| chromVAR | 1.28.0 |
| cicero | 1.3.9 |
| monocle3 | 1.4.26 |

### Genomic and Motif Analysis

| R Package | Version |
|---|---|
| GenomicRanges | 1.58.0 |
| GenomeInfoDb | 1.42.1 |
| IRanges | 2.40.1 |
| rtracklayer | 1.66.0 |
| ChIPpeakAnno | 3.40.0 |
| Biostrings | 2.74.1 |
| BSgenome | 1.74.0 |
| BSgenome.Oryza.sativa.Ensembl.Plants.Osativa.IRGSP | 1.0.0 |
| universalmotif | 1.24.2 |
| TFBSTools | 1.44.0 |
| JASPAR2026 | Not detected |
| RSQLite | 3.53.1 |

### Functional Enrichment

| R Package | Version |
|---|---|
| clusterProfiler | 4.14.6 |
| org.Osativa.eg.db | 1.0.0 |
| enrichplot | 1.26.6 |

### Data Processing

| R Package | Version |
|---|---|
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| tibble | 3.3.1 |
| tidyverse | 2.0.0 |
| reshape2 | 1.4.5 |
| data.table | 1.18.4 |
| Matrix | 1.7-1 |
| purrr | 1.2.2 |
| readxl | 1.5.0 |
| openxlsx | 4.2.8.1 |
| HGNChelper | 0.8.15 |
| BiocParallel | 1.40.0 |
| igraph | 2.3.2 |

### Data Visualization

| R Package | Version |
|---|---|
| ggplot2 | 4.0.3 |
| ggrepel | 0.9.8 |
| patchwork | 1.3.2 |
| pheatmap | 1.0.13 |
| ComplexHeatmap | 2.22.0 |
| circlize | 0.4.18 |
| ggalluvial | 0.12.6 |
| gtable | 0.3.6 |

## 2. Python and Command-Line Software

| Software | Version | Application |
|---|---|---|
| scGLUE | 0.2.3 | Single-cell multi-omics integration |
| deepTools | 3.5.6 | Genomic signal processing and visualization |
| MACS3 | 3.0.2 | Histone modification peak calling |
| Conda | 24.11.3 | Software environment management |

## 3. Reproducibility Notes

The listed versions were obtained from the computational environments used for the Droplet NB-CT analyses. R package versions and Conda-managed software versions may originate from different operating systems or environments.

`JASPAR2026` was not detected in the reported R environment and may need to be installed separately to reproduce analyses that depend on it.

For full reproducibility, users should refer to the analysis scripts and ensure that the appropriate software dependencies, reference genomes, and annotation files are available.

## 4. Related Repository

[Droplet NB-CT GitHub Repository](https://github.com/Tc940417/Droplet-NB-CT)
## Analysis Workflow

### 1. Sequencing Data Processing

Raw DNBC4 sequencing data are processed using `dnbc4tools`. Histone modification peaks are identified using MACS3, followed by construction of peak-by-cell matrices and quality-control analysis with Signac and Seurat.

### 2. Quality Control and Benchmarking (Figure 1)

Scripts in `fig1/` evaluate sequencing quality and data performance, including unique fragments, FRiP, TSS enrichment, sequencing saturation, peak overlap, and correlations with bulk ChIP-seq profiles.

### 3. Cell-Type Identification (Figure 2)

Scripts in `fig2/` and `sctype/` support cell-type annotation using marker genes and single-cell RNA sequencing reference data.

The `GLUE/` directory contains scGLUE analyses integrating snRNA-seq with single-cell H3K4me3 or H3K27ac profiles to facilitate cell-type identification.

### 4. Cell-Type-Specific Chromatin Analysis (Figure 3)

Scripts in `fig3/` characterize cell-type-associated histone modification landscapes through gene activity analysis, peak breadth analysis, differential peak signal visualization, and chromVAR-based transcription factor motif activity analysis.

### 5. Heat-Stress Response Analysis (Figure 4)

Scripts in `fig4/` investigate cell-type-specific responses to heat stress in rice, including:

- Differential H3K4me3 and H3K27ac peaks.
- Differential transcription factor motif activity.
- Integration of histone modification changes with RNA expression.
- Peak-to-gene association and GO enrichment analysis.

## Software Requirements

The main tools and packages used in this repository include:

| Software | Application |
|---|---|
| DNBC4 tools | Sequencing data preprocessing |
| MACS3 | Histone modification peak calling |
| Seurat / Signac | Single-cell data processing and analysis |
| scGLUE | Multi-omics integration |
| chromVAR | Transcription factor motif activity analysis |
| deepTools | Genomic signal analysis |
| scType | Cell-type annotation |
| DESeq2 / MAST | Differential analysis |
| ComplexHeatmap / ggplot2 | Data visualization |
| Cicero | Chromatin co-accessibility analysis |

Additional R and Python packages are specified in the corresponding scripts.

## Data Availability

Raw sequencing data are deposited in the NCBI Sequence Read Archive (SRA).

**BioProject accession:** [PRJNA1476849](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1476849)

Processed datasets are made available separately through GigaDB.

## Citation

If you use the code or analytical workflows from this repository, please cite:

Tong C, et al.

*Droplet NB-CT: A streamlined and highly versatile single-cell chromatin immunoprecipitation method for plants.*

## License

This repository is licensed under the GNU General Public License v3.0 (GPL-3.0). See [LICENSE](LICENSE) for details.
