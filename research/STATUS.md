# Research status

*Last updated: 2026-06-10. This is the one page Spencer needs to read. Everything here
links to a reproducible experiment or a committed document.*

**Goal:** QCE 2027 contributed paper, *"When and why does the homogeneous proxy work?
QAOA parameter setting as subspace compression."* Complete draft by **Aug 15, 2026**;
arXiv late August; submit April 2027. Full plan: see `research/decisions.md` and the
program plan (Claude's plan file, to be mirrored into `research/program.md`).

## Where we are

**Phase 0 (Jun 10–17): formalization + sanity checks.** Just started.

## What we know (established results only)

*(empty — fresh start as of 2026-06-10; results appear here only once they have a
committed, reproducible experiment under `research/experiments/`)*

## Working hypotheses (NOT established — from the deleted research log or intuition)

- H1: The proxy's usefulness is governed by leakage out of the cost-class subspace,
  which concentrates only for ER-like graphs. (The paper's central claim.)
- H2: Fitted proxy shapes (triangle/Gaussian) with lower entrywise MSE against the
  empirical N(c';d,c) gave *worse* parameter setting. (Old log; needs re-verification.)
- H3: Sampling 5–10 bitstrings per cost class suffices to estimate the quantities that
  matter. (Old log; needs re-verification in the leakage metric, not entrywise N.)
- H4: On ER(0.5), random balanced partitions reach ~75–85% approximation ratio, so raw
  AR overstates proxy value. (Internship report; folded into metric design as
  "value-added" = AR_proxy − AR_baseline.)

## Open questions / next experiments

- E0.1 (this week): is the proxy iteration *numerically identical* to the compressed
  statevector evolution (Theorem 1), to 1e-10? Must pass before anything else.
- E1.1–E1.3 (by ~Jul 8): the go/no-go gate — do leakage, fidelity, and parameter-setting
  regret rank graph families identically?

## Decisions needed from Spencer

*(none right now)*
