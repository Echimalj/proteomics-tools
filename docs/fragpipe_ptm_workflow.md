# FragPipe PTM Workflow

This document describes a general workflow for PTM discovery and targeted PTM analysis.

## Overview

The PTM workflow includes:

1. Convert RAW files to mzML
2. Run open or targeted FragPipe searches
3. Filter PSMs for target protein and confidence
4. Summarize modified peptides across conditions
5. Analyze site-level PTM tables
6. Visualize PTM positions along the protein sequence

## Common inputs

```text
psm.tsv
combined_site_STY_79.9663.tsv
combined_site_Kn_42.0106.tsv
```

## Filtering strategy

Recommended filters:

```
target gene/protein only
Q-value <= 0.01
remove decoys
remove contaminants
keep modified peptides only
PTM classification
```

PTMs may be classified as:

Condition A-specific
Condition B-specific
Shared
Notes

Open searches are useful for discovery, but targeted searches are recommended for specific PTMs such as:

phosphorylation
acetylation
ubiquitination
