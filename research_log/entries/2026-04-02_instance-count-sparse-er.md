# Instance Count Sensitivity on Sparse ER

**Date:** 2026-04-02
**Script:** `julia/paper_figures/instance_count_sparse_er.jl`
**Status:** Complete

## Question

How many graph instances does SampN+EmpP need on sparse ER? Sparse ER has
fewer edges per instance, so each provides less information.

## Setup

- n=18, p_edge ∈ {0.2, 0.3, 0.5}
- Instance counts: 1, 3, 5, 10, 20, 40
- 5 random seeds per count, p=1 and p=3

## Key Results

### p=1: 5 instances sufficient

| Instances | ER(0.2) | ER(0.3) | ER(0.5) |
|-----------|---------|---------|---------|
| 1 | 0.736±0.012 | 0.737±0.010 | 0.806±0.008 |
| 5 | 0.775±0.005 | 0.792±0.004 | 0.829±0.005 |
| 10 | 0.774±0.003 | 0.789±0.006 | 0.833±0.005 |
| 20 | 0.776±0.000 | 0.792±0.002 | 0.836±0.000 |

### p=3: 10 instances recommended for sparse ER

| Instances | ER(0.2) | ER(0.3) | ER(0.5) |
|-----------|---------|---------|---------|
| 1 | 0.786±0.012 | 0.792±0.026 | 0.828±0.017 |
| 5 | 0.796±0.029 | 0.842±0.008 | 0.869±0.014 |
| 10 | 0.807±0.010 | 0.830±0.006 | 0.865±0.013 |
| 20 | 0.801±0.009 | 0.824±0.004 | 0.873±0.009 |

## Conclusions

1. **5 instances is sufficient at p=1** for all p_edge values.
2. **10 instances recommended at p=3**, especially for sparse ER (p=0.2)
   where the variance is higher.
3. **20 instances provides minimal additional benefit** over 10.
4. **Sparse ER does need slightly more instances** than ER(0.5), but not
   dramatically more. 10-20 is a safe recommendation.
