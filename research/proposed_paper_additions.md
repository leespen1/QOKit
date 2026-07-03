# Proposed paper additions from the overnight results (E014/E015/E016)

*Written 2026-06-13. These are **ready-to-paste** LaTeX snippets for the three
new results. I did NOT edit §4/§5.4/§7 of `qce2027_paper.tex` directly — the
theory and framing are yours. Accept/reject/edit as you like. Each snippet says
where it goes and which figure to copy into the paper repo's
`Figures/generated/`. All claims trace to a committed experiment.*

---

## 1. E014 → §5.4 (or end of §5.3): argmax transfer, confirmed at depth

This gives §5.4's thesis ("parameter setting is an argmax-transfer problem, not a
state-approximation problem") direct, instance-level evidence at depth, and
explains the "analytical model beats exact compression" surprise.

**Copy:** `research/experiments/014_argmax-robustness/argmax_vs_fidelity.png`
→ `Figures/generated/argmax_vs_fidelity.png`.

```latex
% E014 — argmax transfer is sharpest at depth
The argmax-transfer reading is sharpest at depth. Recomputing the true-landscape
geometry across the seven families ($n=12,14$; $p=1,3$), parameter-setting regret
tracks the proxy's \emph{argmax displacement} from the true optimum at both depths
(Spearman $\rho=0.76$ at $p=1$, $0.63$ at $p=3$), but the state-fidelity deficit
$1-\abs{\braket{\psi}{\phi}}^2$ only at $p=1$: at $p=3$ the proxy-chosen schedules
push fidelity into a narrow band across families and its correlation with regret
collapses to $\rho=-0.02$ (Fig.~\ref{fig:argmax}). Landscape flat-peak robustness
does not predict regret either ($\rho\approx0.1$). The leakage calculus thus bounds
the \emph{state} error, while parameter-setting quality is a \emph{parameter-space}
matter --- which is why a lower-fidelity model whose argmax is better placed (the
analytical $N$ outside dense ER) beats the exact compression. % E014
\begin{figure}[t]
  \centering
  \includegraphics[width=\linewidth]{Figures/generated/argmax_vs_fidelity.png}
  \caption{Parameter-setting regret vs.\ (a) state-fidelity deficit and (b) argmax
  displacement, over $140$ instances at $p=1$ (blue) and $p=3$ (red). Fidelity
  predicts regret at $p=1$ ($\rho=0.39$) but decouples at $p=3$ ($\rho=-0.02$);
  argmax displacement predicts it at both depths ($\rho=0.76,\,0.63$). Parameter
  setting lives in parameter space.} % E014
  \label{fig:argmax}
\end{figure}
```

---

## 2. E015 → §4: the V₂ lemmas (exact, machine-verified)

The ready LaTeX (Lemma "variance reduction" + Lemma "codegree form" + remark) is
in **`research/v2_density_law.md`** (last section). Both lemmas are exact and
verified to 1e-10 on 280 instances; the cubic law $\lambda_1=(|\beta|\gamma^2/8)
\sqrt{V_2}$ holds to <0.3%. They turn §4's "$V_2$ is a within-class variance" into
"$V_2$ is the within-class variance of $\sum_i(\sum_{j\sim i}s_j)^2$, whose
unconditional value is exactly $2\sum_{j\neq k}A_{jk}^2$ (squared codegrees)," and
the Discussion's open $V_2$ problem shrinks to one quantity (the triangle-driven
conditioning correction, confirmed $\propto\tau^2/m$ at Pearson 0.994).

Optional figure: `research/experiments/015_v2-density-law/v2_anatomy.png`
(cubic law on the diagonal + the density-law flattening) → `Figures/generated/`.

---

## 3. E016 → §7 Discussion: the cheap-frame chicken-and-egg does not break

The Discussion currently says instance-adapted low-rank frames "could in principle
beat the homogeneous compression at equal cost, if the chicken-and-egg of needing
states to build the frame can be broken." E016 tests the cheap form and finds it
does not — replace/augment that sentence with:

```latex
% E016 — the chicken-and-egg does not break cheaply
We tested the cheap form of this directly. The depth-$20$ ramp trajectory is
genuinely $\sim\!4$-dimensional (an optimal rank-$4$ frame captures $\approx0.99$
of it), so a good instance-adapted frame would beat the $(m{+}1)$-dimensional
cost-class frame handily. But a rank-$4$ frame built from only the first five
layers is nearly orthogonal to the trajectory's true late-time rank-$4$ subspace
(largest principal angle $71^\circ$--$88^\circ$, growing with ramp) and captures
\emph{less} of the trajectory than the zero-cost cost-class frame. The low-rank
subspace \emph{rotates} across depth, so it cannot be discovered from a cheap
prefix; an instance-adapted frame that beats the cost-class compression would have
to model that rotation explicitly. % E016
```

Optional figure: `research/experiments/016_cheap-prefix-frame/prefix_frame.png`
→ `Figures/generated/`.

---

## Notes

- All three figures are CairoMakie PNGs already committed in their experiment
  directories; copying into `Figures/generated/` is the only step before the
  `\includegraphics` lines resolve.
- E014/E016 are at $n\le14$ (CPU); if you want $n=16$–$20$ versions for the
  paper, they need the HPC GPU path — say the word and I'll write the Slurm driver.
- Update the experiment range to E001--E016 if any of these go in (it currently
  reads E001--E013 after tonight's fix).
