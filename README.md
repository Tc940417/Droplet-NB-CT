# Droplet NB-CT

## Overview

This repository contains the source code used for the analyses presented in the manuscript:

**Droplet NB-CT: A streamlined and highly versatile single-cell chromatin immunoprecipitation method for plants**

The repository provides analysis scripts for single-cell chromatin immunoprecipitation sequencing (NB-CT), chromatin accessibility analysis, multi-omics integration, motif enrichment analysis, and figure generation.

---

## Repository Structure

```
.
├── GLUE/
├── fig1/
├── fig2/
├── fig3/
├── fig4/
├── sctype/
├── Ara-scglue
├── Rice-scglue-H3K27ac
├── Rice-scglue-H3K4me3
├── bulk-pre
├── deeptools
├── macs3
└── LICENSE
```

### Directory description

| Directory | Description |
|------------|-------------|
| GLUE | Scripts for GLUE integration analyses |
| fig1 | Analysis scripts used to generate Figure 1 |
| fig2 | Analysis scripts used to generate Figure 2 |
| fig3 | Analysis scripts used to generate Figure 3 |
| fig4 | Analysis scripts used to generate Figure 4 |
| sctype | Cell type annotation |
| Ara-scglue | Arabidopsis RNA/ATAC integration |
| Rice-scglue-H3K27ac | RNA and H3K27ac integration |
| Rice-scglue-H3K4me3 | RNA and H3K4me3 integration |
| bulk-pre | Bulk sequencing preprocessing |
| deeptools | DeepTools analysis scripts |
| macs3 | Peak calling using MACS3 |

---

## Data Availability

Raw sequencing data have been deposited in the NCBI Sequence Read Archive (SRA).

**BioProject accession:** PRJNA1476849

---

## Software

Major software used in this repository includes:

- R
- MACS3
- deepTools
- SCGLUE
- SCType

Please refer to individual scripts for package dependencies and software versions.

---

## Usage

Each directory contains scripts corresponding to a specific analysis module or figure in the manuscript.

Users may execute the scripts according to the analysis workflow described in the manuscript.

---

## License

All source code in this repository is licensed under the GNU General Public License v3.0 (GPL-3.0).

The directories **fig1–fig4 contain analysis scripts used to generate the corresponding figures and therefore contain source code rather than processed datasets.**

Processed datasets are released separately through GigaDB.

---

## Citation

If you use this repository, please cite:

 Tong C, et al.

Droplet NB-CT: A streamlined and highly versatile single-cell chromatin immunoprecipitation method for plants.
