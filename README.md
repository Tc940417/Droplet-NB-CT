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

## Usage

Clone the repository:

```bash
git clone https://github.com/Tc940417/Droplet-NB-CT.git
cd Droplet-NB-CT
```

Each directory contains scripts associated with a specific computational analysis or figure.

Users should prepare the required input files and modify file paths where necessary before running individual scripts. Some analyses require intermediate Seurat objects, peak matrices, or AnnData files generated in previous steps.

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
