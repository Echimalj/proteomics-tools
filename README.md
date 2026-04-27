# proteomics-tools

Reusable, modular workflows for quantitative proteomics, immunoprecipitation proteomics, and PTM analysis.

🚧 Work in progress

![R](https://img.shields.io/badge/R-4.x-blue)
![FragPipe](https://img.shields.io/badge/FragPipe-MSFragger-orange)
![Status](https://img.shields.io/badge/status-in%20progress-orange)

---

## Overview

This repository provides reusable R utilities for downstream proteomics analysis from FragPipe/MSFragger outputs.

It is designed for workflows involving:

- TMT protein abundance analysis
- Immunoprecipitation enrichment analysis
- IP vs IgG background correction
- genotype- or condition-dependent interactome comparisons
- limma-based statistical modeling
- missing-value handling
- candidate interactor prioritization
- PTM filtering and site-level visualization

The repository is general-purpose. BSN immunoprecipitation proteomics is provided as the main example workflow.

---

## Key Features

**Modular design**  
Reusable functions are separated from project-specific analysis scripts.

**FragPipe-compatible**  
Designed around common FragPipe outputs such as protein abundance tables, PSM tables, and combined PTM site tables.

**Flexible experimental design**  
Supports comparisons such as:

```text
(IP - IgG) in condition A vs (IP - IgG) in condition B
```
**Missing-value aware**
Includes filtering and imputation utilities for proteomics datasets.

**PTM-ready**
Includes utilities for filtering modified peptides and visualizing PTM sites along a protein sequence.

*** Example Use Case: BSN Proteomics
The example scripts use Bassoon (BSN) immunoprecipitation proteomics to demonstrate the workflow.

The main biological comparison is:
(WT_BSN_IP - WT_IgG) - (Mutant_BSN_IP - Mutant_IgG)
This asks which proteins are differentially enriched by BSN pulldown between WT and mutant conditions.

-------
## Repository structure
```
**proteomics-analysis-tools/
├── R/
│   ├── load_fragpipe_utils.R        # Read FragPipe protein/peptide/PTM tables
│   ├── sample_metadata_utils.R      # Build sample metadata and contrasts
│   ├── filtering_utils.R            # Completeness filters and QC summaries
│   ├── imputation_utils.R           # Missing-value imputation helpers
│   ├── limma_utils.R                # limma models and contrasts
│   ├── interactome_utils.R          # IP/IgG enrichment and hit classification
│   ├── ptm_filter_utils.R           # PTM candidate filtering
│   ├── ptm_plot_utils.R             # PTM site/lollipop plotting
│   ├── plotting_utils.R             # Volcano and enrichment plots
│   └── export_utils.R               # Save results and figures
│
├── scripts/
│   ├── bsn_tmt_ip_analysis.R        # Example: BSN IP TMT analysis
│   ├── bsn_interactome_rowmin.R     # Example: BSN row-min interactome workflow
│   └── bsn_ptm_analysis.R           # Example: BSN PTM workflow
│
├── docs/
│   ├── fragpipe_tmt_workflow.md     # General FragPipe TMT workflow notes
│   ├── fragpipe_ptm_workflow.md     # General FragPipe PTM workflow notes
│   └── analysis_notes.md            # Assumptions, thresholds, interpretation
│
├── README.md
└── .gitignore
```

