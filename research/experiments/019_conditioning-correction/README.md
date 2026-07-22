# 019 — Closing the density law: a quadratic bound on the conditioning correction

**Question:** E015 reduced the open V2 problem (the mechanism behind the
empirical density law sqrt(V2) proportional to m) to one quantity: the
conditioning correction Var(E[T|c]) = Var(T) - V2. The crude linear estimate
36 tau^2/m captures the direction but undershoots (measured constant ~60 vs
36). Can the correction be pinned rigorously?

**Answer: Yes, to within a few percent. Projecting T onto span{S, S^2} of
the edge-spin sum S = m - 2c gives a rigorous, polynomial-time lower bound
on Var(E[T|c]) that captures 91-100% (mean 97%) of the exact correction on
all 140 instances across the seven families. Zero bound violations. The
missing channel was quadratic: Cov(T, S^2) = Var(T) exactly, present even
in triangle-free graphs.**

## The result

Four exact moment identities (all verified to 1e-9 per instance, and each a
two-line expansion given E015's codegree lemma T - 2m = 2 sum A_jk s_j s_k):

- Cov(T, S) = 6 tau (each triangle feeds each of its three edges' codegree)
- **Cov(T, S^2) = Var(T) = 2 sum_{j!=k} A_jk^2** (ordered two-paths; the
  surprise that makes the bound tight)
- E[S^3] = 6 tau
- Var(S^2) = 2m^2 - 2m + 24 c4 (c4 = number of 4-cycles)

Since E[T|S] is the L2 projection of T onto all functions of S (and
conditioning on c equals conditioning on S), the projection onto the
subspace span{S, S^2 - m} is a rigorous lower bound:

  Var(E[T|c]) >= v' G^{-1} v,
  v = (6 tau, Var(T)),  G = [[m, 6 tau], [6 tau, 2m^2 - 2m + 24 c4]],

equivalently **V2 <= Var(T) - v' G^{-1} v**, everything computable from the
graph in polynomial time (codegrees, triangles, 4-cycles).

## Measured capture (results.csv, 140 instances, n = 12, 14)

| family | linear only | quadratic bound | V2/Var(T) exact | bound |
|---|---|---|---|---|
| ER(0.5) | 0.64 | 0.96 | 0.352 | 0.376 |
| ER(0.25) | 0.18 | 0.98 | 0.597 | 0.605 |
| BA(k=2) | 0.31 | 0.98 | 0.534 | 0.543 |
| BA(k=4) | 0.61 | 0.96 | 0.356 | 0.382 |
| WS(0.1) | 0.63 | 0.97 | 0.489 | 0.505 |
| WS(0.5) | 0.41 | 0.98 | 0.552 | 0.559 |
| 3-regular | 0.23 | 0.95 | 0.687 | 0.704 |

Pooled capture: mean 0.969, min 0.905, max 1.000; no instance exceeds 1
(consistent with the bound's rigor).

## The ER asymptotic

Substituting ER expectations (m ~ pn^2/2, tau ~ n^3 p^3/6, c4 ~ n^4 p^4/8,
Var(T) ~ 2 n^4 p^4):

  V2 / Var(T)  ->  1 - p (1 + 4p - 2p^2) / (1 + 6p^2 - 4p^3).

At p = 1/2 this gives 0.375 (measured 0.352 at n = 12-14); at p = 1/4 it
gives 0.643 (measured 0.597, larger finite-n gap). The empirical density law
is now an explicit formula up to the <= 9% residual (higher-order channels
S^3, S^4, ...).

## Paper impact

- New Proposition (quadratic conditioning bound) added to section 4 after
  the conditioning remark, with the capture numbers and the ER asymptotic.
- Discussion's open problem updated: what remains is the small residual and
  a rigorous self-averaging statement, not the correction itself.
- Closes journal_readiness gap 5 (density law) in substance.

Reproduce: `JULIA_NUM_THREADS=auto julia --project research/experiments/019_conditioning-correction/run.jl` (~2 min).
