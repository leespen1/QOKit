# Research Log

Structured log of experimental results for the IEEE Quantum Week paper on QAOA
parameter-setting heuristics for non-Erdos-Renyi graphs.

## Protocol

1. **Before running an experiment:** Read `index.md` to check if a similar result
   already exists. Read `next_steps.md` for the highest-priority open question.
2. **After completing an experiment:**
   - Create a new entry in `entries/YYYY-MM-DD_kebab-case-topic.md`
   - Append a row to `index.md`
   - Update `next_steps.md`: mark completed item DONE, add new questions arising

## Operating Modes

**Active mode** (default): Read `next_steps.md`, propose the top P0 item to the
user, wait for approval before executing.

**Autonomous mode** (user says "autonomous" or "overnight"): Pick top P0 item,
execute, log, move to next. When P0 is empty, proceed to P1 items, then P2.
On errors or moderate ambiguity, make a reasonable assumption and proceed, or
skip and move to the next item. Log any skipped items or assumptions in the
entry. Only stop when all priority levels are empty or there is a
critical/unrecoverable issue.

## Entry Template

```markdown
# [Short descriptive title]

**Date:** YYYY-MM-DD
**Paper section:** [see controlled vocabulary below]
**Tags:** [free-form, e.g.: triangle-proxy, erdos-renyi, barabasi-albert, linear-ramp, depth-scaling]
**Status:** [complete | in-progress | superseded-by:FILENAME]

## Motivation
Why this experiment was run. What question it answers. Link to the next_steps.md
item that prompted it, if applicable.

## Setup
- **Script:** `path/to/script.jl`
- **Parameters:** key=value, key=value, ...
- **Output files:**
  - Data: `path/to/data.csv`
  - Figures: `path/to/figure.png`

## Key Findings
1. Most important result first
2. Second finding
3. Third finding

## Significance
How these findings relate to the paper's argument. What they confirm, refute, or
refine about our understanding.

## Next Steps Arising
- [ ] Specific actionable question arising from this result
- [ ] Another question
```

## Controlled Vocabulary: Paper Sections

| Tag | Covers |
|-----|--------|
| `methodology` | Proxy algorithm, homogeneity assumption, proxy fitting |
| `distribution-analysis` | N(c';d,c) shape, Pearson correlations, stddev heatmaps |
| `parameter-transfer` | Transferring parameters across graph sizes/types |
| `approximation-ratios` | QAOA performance comparisons (proxy vs real, proxy types) |
| `scalability` | Depth scaling, large-n behavior, GPU performance |
| `discussion` | Cross-cutting observations, limitations, future work |

## Conventions

- **File naming:** `YYYY-MM-DD_kebab-case-topic.md`. Append `-2`, `-3` for
  multiple entries on the same day/topic.
- **Never delete entries.** Mark superseded entries with `superseded-by:FILENAME`.
- **Index rows** should be under ~120 characters in the summary column.
