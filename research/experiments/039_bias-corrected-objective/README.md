# 023 — A norm-loss bias correction does not help (negative result)

**Question:** E038 showed the normalized proxy overpredicts the objective
everywhere, most at its own argmax. Is that bias predictable from the one
proxy-side observable Theorem 2 gives for free, the norm loss
1 - ||phi||^2, so that a corrected objective F-hat + b (1 - ||phi||^2)
improves parameter setting at zero cost?

**Answer: No. The bias is only weakly and inconsistently correlated with
norm loss (per-instance Pearson -0.26 to 0.60 across the 14 probe grids,
pooled slope +0.053), and applying the pooled linear correction makes
regret WORSE on every instance where it changes anything: p=1 holdout mean
regret 0.047 -> 0.063 (better on 0/126, worse on 112), p=3 holdout 0.043 ->
0.063 (better on 0/126, worse on 87).**

## Reading

Normalization (E036) already removes the norm-visible part of the leakage
bias; what remains is not a function of norm loss, and tilting the
objective by norm loss again just reintroduces a bias with the opposite
sign of the one normalization fixed. This is the third independent
confirmation of the paper's central methodological lesson (after E012's
fitted shapes and E005/E006's filters): calibrating the proxy's VALUES,
by any scalar signal we have tried, does not improve and usually harms its
ARGMAX. The winner's-curse component of the overprediction (E038) is a
selection effect, invisible to any pointwise correction by construction.

## Method

Stage A (probe): instances inst=1 of each (family, n) cell, full p=1 grid;
regress e(theta) = F - F-hat on x(theta) = 1 - ||phi||^2 pooled over the 14
grids (slope 0.0531; intercepts do not move argmaxes). Stage B: corrected
argmax on all 140 instances x p in {1,3}; headline numbers exclude the 14
probe instances.

Reproduce: `JULIA_NUM_THREADS=auto julia --project research/experiments/039_bias-corrected-objective/run.jl` (~12 min).
