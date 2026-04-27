# Proteomics Analysis Notes

## Missing values

Missing values are common in proteomics datasets. This repository includes utilities for:

- completeness filtering
- row-min imputation
- enrichment-based candidate prioritization

Row-min imputation should be interpreted as an exploratory prioritization strategy.

## Interaction contrasts

For IP proteomics, direct IP-vs-IP comparisons can be misleading because they do not account for background binding.

Preferred contrast:

```text
(ConditionA_IP - ConditionA_IgG) - (ConditionB_IP - ConditionB_IgG)
```

## Candidate prioritization

Candidate interactors can be prioritized using:

- effect size
- enrichment above IgG background
- genotype/condition-specific enrichment
- adjusted p-value when sufficiently powered
- PTM interpretation

PTM calls should be interpreted using:

- localization probability
- number of PSMs
- number of runs
- intensity consistency
- biological plausibility
