# 032 — AUDIT: the heavy-tail "proxy wins" niche is a baseline artifact

**Question.** E022–E026 built a positive result — the binned proxy beats
rescaled angle transfer under high-dispersion (Pareto) weights — on a
transfer baseline that rescaled its γ-grid by each instance's **mean edge
weight**. Under heavy tails the mean is tail-dominated and a poor scale.
Does the niche survive a transfer baseline given a *robust* scale?

**Answer. No — the niche is an artifact of the mean-rescaling handicap.**
Headline cell (Pareto(1.5) ER(0.5), n=16, p=1; E022's seeds), mean regret:

| transfer scaled by | mean regret | best-of-3 |
|---|---|---|
| binned proxy | 0.012 | — |
| mean weight (E022) | 0.090 | 0.055 |
| median weight | 0.067 | 0.037 |
| cost-distribution std (fair) | 0.015 | 0.002 |

Rescaling γ by the cost-distribution width — the natural QAOA phase scale,
and what the scale-free quantile-binned proxy implicitly uses — makes
single-source transfer **tie** the proxy and best-of-3 transfer **beat** it.
The proxy does not win a fair fight under heavy tails; it only beat a
transfer baseline handicapped by a scale estimator known to fail exactly in
the tested regime. Consequently the CV²≈1 "boundary" (E023/E026) is
near-tautological: it marks where the mean stops being a usable scale, not
an intrinsic QAOA crossover.

**Scope of the retraction.** This invalidates the *positive* weighted-niche
claim of E022–E026 and its pre-registered-boundary framing. It does NOT
touch: Theorems 1–3 and the error calculus; the V₂ derivation; the
unweighted transfer-beats-proxy results (E018/E019, fair — no scale issue);
E020's binned-compression *machinery* (correct; only the downstream
comparison was unfair); or Max-3-XOR replication.

## Reproduce
```bash
JULIA_NUM_THREADS=8 julia --project research/experiments/032_transfer-scale-audit/run.jl
```
