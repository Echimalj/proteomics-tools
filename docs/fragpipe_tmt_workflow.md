# FragPipe TMT Workflow

This document describes a general workflow for processing TMT proteomics data using FragPipe/MSFragger.

## Overview

The workflow includes:

1. Convert RAW files to mzML using ProteoWizard
2. Search spectra using FragPipe/MSFragger
3. Quantify TMT reporter intensities
4. Export protein-level abundance tables
5. Analyze results in R

## Expected FragPipe output

The main input for downstream R analysis is:

```text
tmt-report/abundance_protein_MD.tsv
```
---

### Example design

A common IP proteomics design is:

```
Condition A IgG
Condition A IP
Condition B IgG
Condition B IP
```
The main contrast is:
```
(ConditionA_IP - ConditionA_IgG) - (ConditionB_IP - ConditionB_IgG)
```
This tests whether IP enrichment differs across conditions.

## Notes 
- Include contaminants in the FASTA database.
- Add decoys in FragPipe.
- For TMT workflows, unused channels should be marked as NA.
- Always inspect known positive controls before downstream analysis.


---
