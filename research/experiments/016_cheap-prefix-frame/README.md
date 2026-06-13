# 016 — Can a cheap instance-adapted low-rank frame beat the cost-class frame?

**Question:** The Discussion calls this "the most interesting open question": the
QAOA ramp trajectory lives in a 2–4-dim *moving* subspace (exp 009), yet the
proxy uses the fixed (m+1)-dim cost-class frame and pays for the mis-aim. Could
an instance-adapted low-rank frame do better — and crucially, can it be built
*cheaply* (from a short prefix), breaking the chicken-and-egg of "needing states
to build the frame"?

**Answer: No, not cheaply — the chicken-and-egg is real. The headroom exists (an
oracle top-4 PCA frame captures ~0.99 of the full p=20 trajectory, far above what
the (m+1)-dim cost-class frame achieves), but a frame built from a cheap 5-layer
prefix is nearly orthogonal to the true late-trajectory subspace (largest
principal angle 71°→83°→88° as the ramp grows toward the working regime) and
captures LESS than the zero-cost cost-class frame in 93–100% of instances (0.70
vs 0.88 at large ramps). The 2–4-dim subspace rotates substantially over depth,
so early cheap states do not reveal the late frame. The cost-class frame, though
provably mis-aimed (exp 009), remains the better zero-knowledge choice; beating
it needs the very states one is trying to avoid computing.**

## Why this matters

Exp 009 showed the cost-class frame is the *wrong* low-dim subspace (mis-aimed,
not too small). The natural follow-up — "then use a better, instance-adapted
frame" — is the paper's headline future-work direction. E016 tests its cheap,
constructive form and finds it does not work: the moving subspace's rotation
across depth defeats prefix-based discovery. This converts the Discussion's
optimism into a concrete, measured obstacle, and sharpens the open problem: a
viable instance-adapted frame must track a *rotating* subspace, which a fixed
cheap-prefix frame cannot.

## Method

Exp 002/004 instance set (same seeds), 7 families × n∈{12,14} × 5 instances ×
3 ramps (small/moderate/large), depth P=20. For each, from the full intermediate
states ψ₀…ψ₂₀ (`qaoa_statevector(...; return_intermediates=true)`):

- **E_cc** — energy captured by the cost-class frame (dim m+1; from
  `project_onto_cost_classes` residuals). Built from costs alone, *zero*
  statevector layers — the proxy's frame.
- **E_full_d** — top-d=4 PCA of the *full* trajectory (oracle ceiling; needs all
  P layers).
- **E_prefix_d** — top-d=4 PCA of only the first K+1=6 layers (the cheap
  constructive frame), then used to capture the *full* trajectory.
- **max principal angle** between the prefix-d and full-d subspaces
  (svdvals(Uᴴ_pre U_full) → cosines; 0° = prefix already spans the late subspace,
  90° = orthogonal).

Capture is `(1/(P+1)) Σ_ℓ ‖Uᴴ ψ_ℓ‖²` for an orthonormal-column frame U.

## Result

Means over instances (n=12,14), depth P=20, d=4, prefix K=5:

| ramp | E_cc (m+1, 0 layers) | E_full_4 (oracle) | E_prefix_4 (cheap) | ∠(prefix,full) | prefix > cost-class |
|---|---|---|---|---|---|
| small    | 0.985 | 0.998 | 0.977 | 71° | 7% of instances |
| moderate | 0.961 | 0.996 | 0.882 | 84° | 0% |
| large    | 0.876 | 0.970 | 0.705 | 88° | 0% |

`prefix_frame.png`: (a) captured energy by the three frames per ramp — the oracle
shows the headroom, the cheap prefix falls below even the cost-class frame at
working ramps; (b) the prefix↔full principal angle, near-orthogonal and growing
with ramp.

**Reading.** The trajectory is genuinely ~4-dim (oracle ≈ 0.99, confirming exp
009), so a good 4-dim instance-adapted frame would beat the cost-class frame
handily. But that subspace is *not* where the early trajectory points: the
prefix frame is 71–88° off and worsens as the ramp (and thus the per-layer
rotation) grows. The mis-aim exp 009 found is therefore not cheaply repairable by
looking at a prefix — the relevant directions only appear late.

## Caveats

- Fixed d=4 and prefix K=5; a longer/adaptive prefix or a larger d would capture
  more (trivially, K→P recovers the oracle), but that defeats "cheap." The point
  is that a *short* prefix at a *small* d does not transfer.
- Capture compares a d=4 frame (prefix/oracle) against the m+1-dim cost-class
  frame — not a fair *dimension* comparison, but the honest *cost* comparison
  (cost-class needs 0 layers; prefix needs K). The decisive, dimension-fair metric
  is the principal angle, which is large regardless.
- n=12,14 only; the subspace-rotation story may sharpen or soften at larger n
  (untested here).
- This is a structural/negative result, not a method. It does not rule out
  cleverer instance-adapted frames (e.g. ones that explicitly model the rotation).

## Reproduce

```
JULIA_NUM_THREADS=auto julia --project research/experiments/016_cheap-prefix-frame/run.jl
julia --project research/experiments/016_cheap-prefix-frame/analyze.jl
julia --project research/experiments/016_cheap-prefix-frame/make_figure.jl
```

Instances seeded by `20260611 + 10000·fam_idx + 100·n + inst` (= exp 002/004/009).
Smoke test: `E16_SMOKE=1 julia --project .../run.jl`.
