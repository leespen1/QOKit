# A plain-language companion to the paper

*Written for Spencer, 2026-06-13. Goal: explain the paper's main ideas in
order, with intuition first and symbols second, so the dense draft
(`papers/OverleafPaper/qce2027_paper.tex`) reads easily afterward. Nothing here
is new science — it is the same content as the paper, unpacked. Where the paper
compresses three ideas into one sentence, this expands them back out.*

---

## 0. The whole paper in one breath

> The homogeneous proxy is **not an approximation of QAOA — it is QAOA with a
> "blur" applied after every layer.** Once you see it that way, everything the
> proxy does, and every way it can fail, becomes measurable with one cheap
> quantity (leakage), and a lot of empirical mystery turns into bookkeeping.

That single reframe — *substitution heuristic → exact projection* — is the
spine of the paper. The three theorems formalize it; the experiments measure
how much the blur costs you in practice; and the "model error" discussion is
about what happens when you don't have the exact blur and have to guess it.

---

## 1. The problem, and what was already known

**QAOA parameter setting is expensive.** To pick good angles
$(\gamma_\ell,\beta_\ell)$ you have to evaluate the objective $\langle C\rangle$
many times, and each evaluation is a full quantum state (on hardware) or a
$2^n$-amplitude simulation (classically). The outer optimization loop multiplies
that cost into every step.

**The homogeneous proxy (Sud et al. 2024) is a shortcut.** Its bet: bitstrings
that have the *same cost* tend to pick up *nearly the same amplitude* during
QAOA. If that were exactly true, you wouldn't need $2^n$ amplitudes — you'd need
only **one amplitude per distinct cost value** (at most $m+1$ of them for an
$m$-edge MaxCut graph). The proxy tracks exactly that: a short vector
$Q_\ell(v)$, one complex number per cost $v$, evolved by a recursion
(Eq. 2/8) whose coefficients are the **cost-and-distance distribution**
$N(v';d,v)$ — "starting from a typical bitstring of cost $v'$, how many
bitstrings of cost $v$ sit at Hamming distance $d$?"

**What was left open** (this is the paper's hook):

1. *What is the proxy, mathematically?* Sud et al. derived it by substituting a
   class-averaged $N$ into a sum-over-paths expansion, after which — in their own
   words — the evolution "is no longer restricted to unitary evolution" and the
   $Q_\ell(v)$ are "analogues of amplitudes." That's a **description of how it's
   computed, not a characterization of what it is.**
2. *When does it work?* They validated it only on Erdős–Rényi (ER) MaxCut. Is
   the mechanism ER-specific? What makes it degrade with depth or with bigger
   angles? Nobody knew.

The paper answers #1 *exactly* (a theorem) and #2 *quantitatively* (an error
calculus + experiments).

---

## 2. The key reframe — Theorem 1 (exactness)

Here is the picture to hold in your head.

Take the true QAOA state, a vector of $2^n$ amplitudes. Group the bitstrings
into **cost classes**: bucket $S_v$ holds every bitstring with cost $v$. Now do
this operation, call it $P$:

> **Within each bucket, replace every amplitude by the bucket's average.**

That's it. $P$ "flattens" the state so it's constant on each cost class — it
*blurs away* any variation between bitstrings that share a cost. A state that
survives $P$ unchanged (already constant on every class) is called **Perfectly
Homogeneous**; these states form a subspace $\mathcal{H}_{\mathrm{hom}}$ of
dimension $\le m+1$. $P$ is the orthogonal projector onto it.

**Theorem 1 says:** running the proxy for $\ell$ layers gives *exactly* the same
thing as running true QAOA but inserting this blur $P$ after every single layer:

$$\ket{\phi_\ell} \;=\; P\,U_\ell\,P\,U_{\ell-1}\cdots P\,U_1\,\ket{+}^{\otimes n}.$$

The proxy's mysterious "analogue of amplitude" $Q_\ell(v)$ is literally **the
common amplitude that all bitstrings in bucket $v$ have after you blur.** No
approximation has been made *in the definition* — the proxy step and
"evolve-one-layer-then-blur" are the same arithmetic. (Experiment 001 confirms
this numerically to $2.5\times10^{-16}$, i.e. machine precision.)

**Why this matters — three immediate payoffs:**

- **The non-unitarity is demystified.** The proxy isn't unitary because *a
  projection isn't unitary* — it throws away the part of the state that stuck
  out of $\mathcal{H}_{\mathrm{hom}}$. The "lost norm" is not a bug; it is a
  **measurable error budget** (see Prop. / §3).
- **You get a free error certificate.** Blurring can only *shrink* a vector
  (projections are contractions), so the proxy's norm can only go *down*. If you
  ever run a proxy and its norm goes *up*, the $N$ you used cannot be the real
  instance's $N$. Norm inflation = proof of model error, at zero cost. (This
  becomes important in §5 — the analytical model trips it on dense graphs.)
- **It splits the error cleanly** — which is the next section, and the single
  most important conceptual tool in the paper.

> **One-line version of Theorem 1:** *The proxy = QAOA + blur-after-each-layer.
> The blur is an orthogonal projection onto "constant-on-cost-classes" states.*

---

## 3. The two kinds of error (read this twice)

Theorem 1 holds **only when $N$ is the instance's own exact empirical average**
(Eq. 1). But nobody uses that in practice — computing it is $O(4^n)$, as
expensive as the thing you're trying to avoid. In practice you use an
*approximate* $N$: the analytical binomial/multinomial formula of Sud et al., or
a fitted triangle/Gaussian shape, or a sampled estimate. Each choice replaces
the exact transfer matrix $T$ with some $T_{\text{model}} = T + E$.

So the **total error of a real proxy run = two independent pieces:**

| Error | Where it comes from | Controlled by |
|---|---|---|
| **Compression error** | The blur $P$ itself, even with perfect $N$ — the proxy lives in an $(m{+}1)$-dim subspace and QAOA doesn't. | **Leakage** (§4). Model-independent. |
| **Model error $E$** | Using an approximate $N$ instead of the instance's true one. | How good your distribution model is. |

**The paper's recurring punchline is that these two behave independently** and
must be diagnosed differently. The old G-RIPS worry — "the proxy fails on
non-ER graphs" — is, the paper shows, a *model-error* story (the analytical
formula is wrong off ER), **not** a compression-error story (the blur itself is
near-lossless on every family at $p=1$). Conflating them is exactly the mistake
that made the proxy look mysterious.

Keep this table in mind for the rest of the paper. §3–4 are entirely about the
*left* column (compression error). §5.4 is about the *right* column (model
error).

---

## 4. The error calculus — Theorems 2 and 3

Now that compression error has a name, the paper builds three cheap tools to
measure and predict it.

### 4a. Leakage and the telescoping bound (Theorem 2)

Define the **leakage of layer $\ell$**:

$$\lambda_\ell \;=\; \big\lVert (I-P)\,U_\ell\,\ket{\phi_{\ell-1}} \big\rVert.$$

In words: take the (already-blurred) proxy state, apply one true QAOA layer, and
**measure how much of the result sticks out of the homogeneous subspace** — the
part the next blur will throw away. That's the amplitude that "stopped being a
function of cost" during this layer.

Two facts (Theorem 2), both one-line proofs:

- **Telescoping bound:** $\displaystyle \lVert \psi_p - \phi_p\rVert \le
  \sum_{\ell=1}^p \lambda_\ell.$ The total distance between true QAOA and the
  proxy is **at most the per-layer leakages added up.** Errors accumulate at
  worst linearly.
- **Norm bookkeeping:** $\lVert\phi_\ell\rVert^2 = \lVert\phi_{\ell-1}\rVert^2 -
  \lambda_\ell^2$ (Pythagoras). The norm you lose each layer *is* the leakage
  squared.

**Why this is the workhorse of the paper:** $\lambda_\ell$ costs **one
statevector layer + one $O(2^n)$ averaging pass** to measure. No $O(4^n)$
distribution, no $2^n\times 2^n$ matrices. So leakage is a *practical
instrument* — you can compute it, plot it, and use it to rank graph families.
Experimentally the bound's slack is only $3$–$6\times$ at $p=20$ (never
vacuous), and $\sum\lambda_\ell$ is a near-functional *predictor* of the actual
error, not just an upper bound.

### 4b. What leakage actually is — Theorem 3 (variance identity)

Theorem 3 is the conceptual heart. It says the per-class leakage equals a
**within-class variance**:

$$\lambda(v')^2 \;=\; \frac{1}{M_{v'}}\sum_{v} M_v\,\operatorname{Var}_{y\in
S_v}\!\big[g_{v'}(y)\big], \qquad g_{v'}(y)=\sum_d f_d(\beta)\,n(y;d,v').$$

Unpack it: $g_{v'}(y)$ is the **mixer-weighted neighborhood profile** of
bitstring $y$ — "how strongly does the mixer connect $y$ back to cost-$v'$
states, summed over distances, weighted by $f_d(\beta)$." Leakage is large
exactly when **bitstrings that share a cost $v$ have *different* $g_{v'}$
values** — i.e. when "same cost" fails to imply "same neighborhood." That is the
precise, quantitative version of the proxy's founding intuition ("the proxy is
good when bitstrings of equal cost have equal neighborhoods").

The crucial subtlety — and a genuinely useful result — is the **weighting**.
$\lvert f_d(\beta)\rvert = (\cos\beta)^{n-d}(\sin\beta)^d$ falls off fast in $d$
for small $\beta$, so **only the small-distance (local) part of the neighborhood
variance matters** in the regime where the proxy is actually used.

> **Consequence for the G-RIPS triangle/Gaussian fitting:** fitting a shape to
> $N$ by *unweighted entrywise MSE* optimizes the **wrong norm.** It spends
> effort matching large-$d$ entries that $f_d(\beta)$ multiplies by nearly zero.
> This is the theoretical seed of the "fitted-shape paradox" in §5.4.

### 4c. Why it's so accurate at small angles — the Lemma + cubic Corollary

This is the prettiest part and worth the five minutes.

**Step 1: half of every layer never leaks.** A QAOA layer is
$U = B(\beta)\Phi(\gamma)$. The phase separator $\Phi(\gamma)=e^{-i\gamma C/2}$
is diagonal in cost, so it multiplies each class state $\ket{v}$ by a single
phase $e^{-i\gamma v/2}$ — it maps homogeneous states to homogeneous states
**exactly, on every graph, at every angle.** *All leakage comes from the mixer
$B(\beta)$.* (This is why you'll see the paper say "inhomogeneity is injected
only by the mixer.")

**Step 2: the mixer's first-order leak cancels — and that's a MaxCut identity,
not luck.** Expand the mixer to first order: $B \approx I - i\beta\sum_j X_j$.
The $\sum_j X_j$ part connects $y$ to its $n$ single-bit-flip neighbors. To know
how much this leaks, you need the distribution of *neighbor costs* $c(y\oplus
e_i)$. The **Lemma** computes its first moment exactly:

$$\sum_{i=1}^n c(x\oplus e_i) = (n-4)\,c(x) + 2m \quad\text{for every } x.$$

The sum of all single-flip-neighbor costs is a **fixed linear function of
$c(x)$** — so it is *the same for every bitstring of a given cost*, i.e.
class-constant, i.e. **annihilated by $(I-P)$.** (The proof is one line: flipping
bit $i$ toggles the edges at vertex $i$; sum over $i$ using $\sum_i\deg(i)=2m$.
It's a consequence of MaxCut being **2-local** — every edge has two
endpoints — not of any randomness.)

**Step 3: the survivor is third order.** With the constant ($k{=}0$) and linear
($k{=}1$) moments killed, the first term that leaks is the **second moment** of
neighbor costs. The result (**Corollary**):

$$\lambda_1 = \frac{|\beta|\,\gamma^2}{8}\sqrt{V_2}\,\big(1+O(\gamma)+O(\beta)\big),
\qquad V_2 = 2^{-n}\sum_v M_v\operatorname{Var}_{S_v}[s_2],\; s_2(y)=\textstyle\sum_i c(y\oplus e_i)^2.$$

So compression error is **third order in the angles** ($\beta\gamma^2$), with a
prefactor that is a within-class variance of *squared* neighbor costs. Measured
exponents: $0.98$ in $\beta$, $1.92$ in $\gamma$ — i.e. $1$ and $2$. This is the
quantitative explanation of the long-observed "the proxy is accurate at small
$\gamma$" — now with the mechanism, the constant, and the $\beta$-dependence
spelled out.

**Bonus consequence:** for a linear-ramp schedule with fixed endpoints, each
$\lambda_\ell$ depends only on $(\gamma_\ell,\beta_\ell)\approx f(\ell/p)$, so
$\sum_\ell\lambda_\ell \approx p\cdot\overline{f}$ — total leakage grows
**linearly in depth $p$.** Predicted ratio $30/20 = 1.5$; measured $1.505$.

### 4d. The §4.3 "ladder" — why homogeneity happens at all

§4.3 assembles the above into a narrative answer to "why do equal-cost
bitstrings end up with equal amplitudes?":

1. The initial state $\ket{+}^{\otimes n}$ is **exactly** homogeneous.
2. The phase separator preserves homogeneity **exactly** (Step 1).
3. The mixer is the only source of leakage, and Theorem 3 says only
   *within-class variance* of the weighted neighborhood enters.
4. At small $\beta$ that variance is dominated by single flips, where the Lemma
   makes the first moment an exact function of cost — **leading-order protection
   on any graph** (it's 2-locality, not randomness).
5. Randomness buys only the *next* rung: $V_2$ self-averages over the
   exponentially large cost classes for random-like instances. **This is the one
   rung that is measured, not proved** — and it's exactly the open "$V_2$
   problem" flagged in the discussion.
6. Finally it's partly self-fulfilling: QAOA's good schedules shrink $\gamma$
   with density, so the angles worth running are precisely the low-leakage ones.

---

## 5. What the experiments actually show (organized by question)

The §5 prose is a wall of numbers; here it is reorganized as the questions each
result answers.

**Q: Is compression error about being Erdős–Rényi, or something simpler?**
→ **Density, full stop.** Family-mean leakage tracks edge count at Pearson
$0.97$; $\lambda/(\beta\gamma^2 m)$ is constant to $\pm40\%$ across all seven
families. Surprisingly, ER(0.5) — the family the heuristic was built for — is
the *worst*-compressing; random 3-regular the best. (Caveat: this fixes the
angles; shrinking $\gamma$ with density absorbs part of it.)

**Q: Is the telescoping bound actually usable, or just true?**
→ **Usable.** Slack $3$–$6\times$ at $p=20$ over 840 runs; near-functional
predictor; linear-in-depth ($1.505$ vs $1.5$).

**Q: When the proxy degrades, is it because the state got complicated, or
because the fixed cost-class frame points the wrong way?**
→ **Wrong way (mis-aim), not complexity.** Stacking $p=20$ trajectory states and
doing PCA shows effective dimension $2$–$4$. The state stays simple; the *fixed*
cost-class subspace just drifts off the tiny moving subspace the trajectory
actually lives in. (This is the most interesting open door — see §7.)

**Q: Does low compression error actually mean good parameters?**
→ **At $p=1$, yes, and the ranking holds:** regret $0.03$–$0.06$ at $p=1$,
positive value-added on every family, mild growth with $n$. Ranking families by
leakage reproduces their ranking by regret (Spearman $\rho=0.96$ at $n{=}12$,
$0.86$ at $n{=}14$). **Honest limit:** at deeper $p$ and larger $n$ the regret
differences between families shrink to a few $0.01$ and the ranking signal fades
with them. (This is the paper's central claim and also its most fragile one —
see the critique.)

**Q (model error): how bad is using the analytical $N$ off-ER?**
→ **Its argmax is excellent, its values are garbage.** The analytical formula
(with an effective edge probability) gets regret $0.01$–$0.02$ on
BA/WS/3-regular/sparse-ER — *better* than exact compression, because the smooth
model regularizes instance noise. But its predicted $\langle C\rangle$ is
meaningless in absolute terms, and on **dense ER(0.5)** its argmax lands on an
unphysical spurious peak ($\langle C\rangle=93$ on a 38-edge graph, norm
inflated $7.5\times$). The norm certificate *detects* this but **cannot repair
it** — rejecting norm-inflated predictions everywhere also rejects the (correct)
answer on sparse graphs. Reported as a negative result (experiments 004–006).

**Q: does fitting better shapes to $N$ help?**
→ **No — the fitted-shape paradox.** Fitting a triangle/Gaussian to $N$ by
entrywise MSE *improves the fit* yet *worsens or doesn't change* parameter
setting (triangle regret $0.054\to0.211$, worse on 136/140; Gaussian improves
MSE $10\times$ and changes the chosen angles on $0/140$). **No scalar
mismatch-norm we tested predicts regret** (entrywise MSE $\rho\approx-0.1$,
amplitude error $\rho\approx-0.2$, landscape correlation $\rho\approx+0.1$) —
only **argmax displacement** does ($\rho\approx0.7$). Lesson: parameter setting
is an *argmax-transfer* problem, not a *state-approximation* problem. (This is
Theorem 3's wrong-norm point, confirmed from the empirical side.)

**Q: so what should a practitioner actually do?**
→ **Sampled $N$ ($S\approx10$ bitstrings per cost class) + the empirical cost
distribution.** Matches exact parameter setting within $0.02$ AR, regret no
worse, family-agnostic, no dense-graph artifact. **Honest scope:** its cost is
$O(S\,m\,2^n)$ — still exponential, useful only where statevector passes are
affordable but the $O(4^n)$ exact distribution is not (roughly $n\lesssim30$).
Whether $N$ has a polynomial Monte-Carlo estimator is open; beyond classical
reach, the analytical model is the only option, and there you lean on its argmax
robustness and avoid dense graphs.

---

## 6. What is proved vs. measured (so you can defend it)

- **Proved (exact, machine-verified):** Theorem 1 (proxy = compression),
  Theorem 2 (telescoping bound + norm identity), Theorem 3 (variance identity),
  the Lemma (neighbor-cost sum), the Corollary's cancellation (first-order leak
  is identically zero). These hold for *any* graph; no randomness needed.
- **Measured, not proved:** that $V_2$ self-averages so that leakage tracks
  density across random-graph families; the family rankings; all regret numbers;
  the sampled-$N$ recipe's quality. These are empirical over the tested range
  ($n\le20$, $p\le30$, 7 families, 20–30 instances/cell, fixed seeds).
- **The honest seams** (the paper states all of these):
  - The leakage→regret ranking is a **7-point Spearman** that is only resolvable
    at $p=1$/small $n$; it fades at scale. Don't oversell it.
  - Everything is **unweighted MaxCut**, $n\le20$.
  - Trajectory PCA is a *diagnostic*, not a method (it needs the states the proxy
    exists to avoid computing).
  - The recipe's cost is still exponential.

---

## 7. The two open doors (for "what's next")

1. **Instance-adapted low-rank frames.** §5 shows the trajectory lives in a
   $2$–$4$-dim *moving* subspace, but the proxy uses a *fixed* $(m{+}1)$-dim
   cost-class frame and pays for the mismatch. Could a cheap short-depth state
   give a better frame and beat the homogeneous compression at equal cost? The
   chicken-and-egg ("you need states to build the frame") is the obstacle.
2. **An analytic $V_2$ per family.** Turning the measured density law
   ($\lambda/(\beta\gamma^2 m)\approx$ const) into a *theorem* by computing the
   expected within-class variance of $s_2$ for each random-graph ensemble. The
   discussion calls this "open and tractable." (I'm attempting a first pass
   tonight — see the research log.)

---

## 8. Symbol cheat-sheet (friendlier than the paper's Table I)

| Symbol | Say it as | Meaning |
|---|---|---|
| $n,\,m,\,p$ | "qubits, edges, depth" | vertices; edges (= max cost); QAOA layers |
| $c(x)$ | "cost of $x$" | cut value of bitstring $x$ |
| $S_v,\,M_v$ | "cost class $v$, its size" | all $x$ with $c(x)=v$; how many there are |
| $\ket{v}$ | "class state" | uniform superposition over $S_v$ (normalized) |
| $\mathcal{H}_{\mathrm{hom}}$ | "homogeneous subspace" | states constant on each cost class ($\dim\le m{+}1$) |
| $P$ | "the blur" / "class averaging" | projector onto $\mathcal{H}_{\mathrm{hom}}$: average within each bucket |
| $U_\ell=B(\beta)\Phi(\gamma)$ | "a QAOA layer" | mixer $\times$ phase separator |
| $\Phi(\gamma)$ | "phase separator" | diagonal; **preserves homogeneity exactly** |
| $B(\beta)$ | "mixer" | the **only** source of leakage |
| $n(x;d,v)$ | "neighborhood profile" | # bitstrings at distance $d$ from $x$ with cost $v$ |
| $N(v';d,v)$ | "homogeneous distribution" | class-average of $n$; the recursion's coefficients |
| $Q_\ell(v)$ | "proxy amplitude" | the one number per cost the proxy tracks $=$ blurred amplitude |
| $\lambda_\ell$ | "leakage" | norm that leaves $\mathcal{H}_{\mathrm{hom}}$ in layer $\ell$ |
| $\phi_\ell$ vs $\psi_\ell$ | "compressed vs true" | proxy(=blurred) trajectory vs real QAOA |
| regret | "regret" | (best AR on the grid) − (AR at the proxy's chosen angles) |
| value-added | "value-added" | (proxy AR) − (random balanced-partition AR ≈ 0.75) |
| argmax transfer | "does the peak land right?" | the only thing parameter setting actually consumes |

---

## If you remember five things

1. **Proxy = QAOA + a blur (projection) after every layer** (Thm 1). Not an
   approximation by definition; an exact compression.
2. **Two errors, independent: compression error (the blur) and model error (a
   wrong $N$).** Most past confusion is mixing them up.
3. **Leakage** $\lambda_\ell$ is the cheap ($O(2^n)$) ruler for compression
   error; it bounds the true–proxy distance (Thm 2) and equals a *mixer-weighted
   within-class variance* (Thm 3).
4. **At small angles, compression error is third order ($\beta\gamma^2$)**
   because the phase separator never leaks and the mixer's first-order leak
   cancels by a MaxCut identity. Density (edge count) sets the constant.
5. **Parameter setting is argmax-transfer, not state-approximation.** That's why
   fitting $N$ by MSE backfires, why the analytical model works far outside its
   domain, and why "sampled $N$ + empirical $P$" is the recommended recipe.
