# A plain-language companion to the paper

*(Note: this file is the companion to the **paper** specifically. The living
walkthrough of all results, including experiments not yet folded into the
paper, is [explainer.md](explainer.md).)*

*Written for Spencer, 2026-06-13; rewritten 2026-07-31 against the
reconciled QINP manuscript. The canonical paper is now
`papers/OverleafPaper/qinp_paper.tex` (single-column Springer QINP journal
manuscript, ~20 pp; `qce2027_paper.tex` is the 9 pp conference cut of the
same content). On 2026-07-31 the two diverged work lines were reconciled:
the July editorial line is canonical, and the loop-1 science was ported in
with experiment refs renumbered E033--E043 (matching
`research/experiments/033_*..043_*`); E017--E032 are the July series
(`017_error-directions` .. `032_transfer-scale-audit`); E001--E016 are
shared by both lines. Goal of this file: explain the paper's ideas in
order, intuition first and symbols second, completely enough that you could
rewrite the paper from scratch from this document alone. Nothing here is
new science; it is the paper's content, unpacked. Section 10 gives the
paper's skeleton, section 11 the key-numbers table with experiment
provenance, and section 8 marks the one frontier (E041--E043) that is NOT
yet in the paper.*

---

## 0. The whole paper in one breath

> The homogeneous proxy is **not an approximation of QAOA. It is QAOA with a
> "blur" applied after every layer.** Once you see it that way, everything
> the proxy does, and every way it can fail, becomes measurable with one
> cheap quantity (leakage), and a lot of empirical mystery turns into
> bookkeeping. And when you then calibrate the proxy honestly against the
> cheapest baselines, it turns out **not to be a faster or better parameter
> setter**: plain within-family angle transfer is at least as good as every
> proxy variant tested, at every depth and weight regime given a fair scale.
> The lasting contribution is the **error calculus itself**: a cheap,
> measurable certificate of how much of the QAOA state any cost-class
> surrogate silently discards.

Two reframes, then, form the spine: *substitution heuristic to exact
projection* (the theorems), and *speedup claim to error instrument* (the
experiments plus the honest baseline calibration). The paper's title
question, "when does the homogeneous proxy work?", gets a quantitative
answer: the compression is near-lossless wherever density-controlled leakage
is small; parameter quality is an argmax-transfer matter that decouples from
state fidelity; and against well-scaled transfer the proxy matches but never
beats, extending Sud et al.'s own low-depth parity finding to high depth.

---

## 1. The problem, and what was already known

**QAOA parameter setting is expensive.** To pick good angles
$(\gamma_\ell,\beta_\ell)$ you have to evaluate the objective
$\langle C\rangle$ many times, and each evaluation is a full quantum state
(on hardware) or a $2^n$-amplitude simulation (classically). The outer
optimization loop multiplies that cost into every step.

**The homogeneous proxy (Sud et al. 2024) is a shortcut.** Its bet:
bitstrings that have the *same cost* tend to pick up *nearly the same
amplitude* during QAOA. If that were exactly true, you wouldn't need $2^n$
amplitudes; you'd need only **one amplitude per distinct cost value** (at
most $m+1$ of them for an $m$-edge MaxCut graph). The proxy tracks exactly
that: a short vector $Q_\ell(v)$, one complex number per cost $v$, evolved
by a recursion whose coefficients are the **cost-and-distance distribution**
$N(v';d,v)$: "starting from a typical bitstring of cost $v'$, how many
bitstrings of cost $v$ sit at Hamming distance $d$?"

For the rewrite you need the recursion itself. With mixer matrix elements
$f_d(\beta)=(\cos\beta)^{n-d}(-i\sin\beta)^d$ depending only on Hamming
distance $d$:

$$Q_\ell(v') = \sum_{d,v} f_d(\beta_\ell)\, e^{-i\gamma_\ell v/2}\,
Q_{\ell-1}(v)\, N(v';d,v), \qquad Q_0 \equiv 2^{-n/2},$$

and the empirical $N$ is the class average
$N(v';d,v)=\frac{1}{M_{v'}}\sum_{x\in S_{v'}} n(x;d,v)$. After $p$ layers
the predicted objective is $\langle C\rangle \approx \sum_{v}
M_v\,|Q_p(v)|^2\,v$ (one term per attained cost; the paper never needs the
$2^n\times 2^n$ picture). Each layer costs $O(nm^2)$.

**What Sud et al. actually claimed** (the paper now takes both claims as
its starting point and completes the comparison they left open):

- On $G(20,\tfrac12)$ at $p\le3$, their heuristic was **statistically
  indistinguishable** from median-parameter transfer (mean per-instance
  difference $\le0.02$ AR, mixed in sign). They reported parity, not
  superiority, at low depth.
- Restricting to 4-parameter linear ramps, monotonically increasing AR out
  to $p=20$, a regime for which **no transfer table then existed** (their
  free-parameter transfer loop did not converge at depth).

**What was left open** (this is the paper's hook):

1. *What is the proxy, mathematically?* Sud et al. derived it by
   substituting a class-averaged $N$ into a sum-over-paths expansion, after
   which, in their own words, the evolution "is no longer restricted to
   unitary evolution" and the $Q_\ell(v)$ are "analogues of amplitudes."
   That is a **description of how it's computed, not a characterization of
   what it is.**
2. *When does it work?* They validated it only on Erdos-Renyi (ER) MaxCut.
   Is the mechanism ER-specific? What makes it degrade with depth or with
   bigger angles? Nobody knew. And at $p=10$--$20$, does the heuristic
   still beat transfer once ramp transfer exists?

The paper answers #1 *exactly* (a theorem) and #2 *quantitatively* (an
error calculus, experiments, and the completed head-to-head).

---

## 2. The key reframe: Theorem 1 (exactness)

Here is the picture to hold in your head.

Take the true QAOA state, a vector of $2^n$ amplitudes. Group the bitstrings
into **cost classes**: bucket $S_v$ holds every bitstring with cost $v$. Now
do this operation, call it $P$:

> **Within each bucket, replace every amplitude by the bucket's average.**

That's it. $P$ "flattens" the state so it's constant on each cost class; it
*blurs away* any variation between bitstrings that share a cost. A state
that survives $P$ unchanged (already constant on every class) is called
**Perfectly Homogeneous**; these states form a subspace
$\mathcal{H}_{\mathrm{hom}}$ of dimension $\le m+1$. $P$ is the orthogonal
projector onto it.

**Theorem 1 says:** running the proxy for $\ell$ layers gives *exactly* the
same thing as running true QAOA but inserting this blur $P$ after every
single layer:

$$\ket{\phi_\ell} \;=\; P\,U_\ell\,P\,U_{\ell-1}\cdots P\,U_1\,\ket{+}^{\otimes n}.$$

The proxy's mysterious "analogue of amplitude" $Q_\ell(v)$ is literally
**the common amplitude that all bitstrings in bucket $v$ have after you
blur.** No approximation has been made *in the definition*; the proxy step
and "evolve-one-layer-then-blur" are the same arithmetic. (Proof: a
two-line induction. Group the mixer sum by distance and cost; class
averaging turns $n(x;d,v)$ into $N(v';d,v)$. Experiment 001 confirms it
numerically to $2.5\times10^{-16}$, machine precision. The equivalent
transfer-matrix form is $T = D^{-1/2}(PUP)|_{\mathcal{H}_{hom}}D^{1/2}$
with $D=\mathrm{diag}(M_v)$.)

**Why this matters, three immediate payoffs:**

- **The non-unitarity is demystified.** The proxy isn't unitary because *a
  projection isn't unitary*: it throws away the part of the state that
  stuck out of $\mathcal{H}_{\mathrm{hom}}$. The "lost norm" is not a bug;
  it is a **measurable error budget** (section 4).
- **You get a free error certificate** (the paper's Proposition 1).
  Blurring can only *shrink* a vector (projections are contractions), so
  the proxy's norm can only go *down*. If you ever run a proxy and its norm
  goes *up*, the $N$ you used cannot be the real instance's $N$. Norm
  inflation = proof of model error, at zero cost. (The analytical model
  trips it on dense graphs, section 5.) The paper immediately flags that
  **the certificate's converse fails in an instructive way**: no inflation
  does not mean no model error, and, sharper (E017), rejecting inflated
  predictions or dividing them out can *discard* argmax signal. That
  refinement recurs throughout: the norm machinery is an *instrument*, not
  a recipe.
- **It splits the error cleanly**, which is the next section and the single
  most important conceptual tool in the paper.

Three delimiting remarks the paper attaches to Theorem 1: it holds for the
same-instance empirical $N$ only (any analytical, averaged, or fitted $N$
defines a different recursion $T_{\mathrm{model}}=T+E$); unattained costs
are inert (the empirical $N$ zeroes them after one step); and the recursion
must not renormalize (that is what makes the norm a certificate).

> **One-line version of Theorem 1:** *The proxy = QAOA + blur-after-each-layer.
> The blur is an orthogonal projection onto "constant-on-cost-classes" states.*

---

## 3. The two kinds of error (read this twice)

Theorem 1 holds **only when $N$ is the instance's own exact empirical
average**. But nobody uses that in practice; computing it is $O(4^n)$, as
expensive as the thing you're trying to avoid. In practice you use an
*approximate* $N$: the analytical binomial/multinomial formula of Sud et
al., or a fitted triangle/Gaussian shape, or a sampled estimate. Each
choice replaces the exact transfer matrix $T$ with some
$T_{\text{model}} = T + E$.

So the **total error of a real proxy run = two independent pieces:**

| Error | Where it comes from | Controlled by |
|---|---|---|
| **Compression error** | The blur $P$ itself, even with perfect $N$; the proxy lives in an $(m{+}1)$-dim subspace and QAOA doesn't. | **Leakage** (section 4). Model-independent. |
| **Model error $E$** | Using an approximate $N$ instead of the instance's true one. | How good your distribution model is. |

**The paper's recurring punchline is that these two behave independently**
and must be diagnosed differently. The old G-RIPS worry, "the proxy fails
on non-ER graphs," is, the paper shows, a *model-error* story (the
analytical formula is wrong off ER), **not** a compression-error story (the
blur itself is near-lossless on every family at $p=1$). Conflating them is
exactly the mistake that made the proxy look mysterious.

Keep this table in mind. Section 4 is entirely about the *left* column.
The model-error subsection of the experiments (5.4 in the paper) is about
the *right* column.

---

## 4. The error calculus: Theorems 2 and 3, and the MaxCut structure

Now that compression error has a name, the paper builds cheap tools to
measure and predict it.

### 4a. Leakage and the telescoping bound (Theorem 2)

Define the **leakage of layer $\ell$**:

$$\lambda_\ell \;=\; \big\lVert (I-P)\,U_\ell\,\ket{\phi_{\ell-1}} \big\rVert.$$

In words: take the (already-blurred) proxy state, apply one true QAOA
layer, and **measure how much of the result sticks out of the homogeneous
subspace**, the part the next blur will throw away. That's the amplitude
that "stopped being a function of cost" during this layer.

Two facts (Theorem 2), both one-line proofs:

- **Telescoping bound:** $\displaystyle \lVert \psi_p - \phi_p\rVert \le
  \sum_{\ell=1}^p \lambda_\ell.$ The total distance between true QAOA and
  the proxy is **at most the per-layer leakages added up.** Errors
  accumulate at worst linearly.
- **Norm bookkeeping:** $\lVert\phi_\ell\rVert^2 =
  \lVert\phi_{\ell-1}\rVert^2 - \lambda_\ell^2$ (Pythagoras). The norm you
  lose each layer *is* the leakage squared.

**Why this is the workhorse:** $\lambda_\ell$ costs **one statevector layer
plus one $O(2^n)$ averaging pass** to measure. No $O(4^n)$ distribution, no
$2^n\times 2^n$ matrices. So leakage is a *practical instrument*.
Experimentally the bound's slack at $p=20$ is about $4\times$ in the median
and never above $5.8\times$ over 840 runs, and $\sum\lambda_\ell$ is a
near-functional *predictor* of the actual error, not just an upper bound
(E003). The paper is careful about novelty here: the proof technique is the
standard one for projected dynamics (a discrete-time, fixed-subspace
analogue of the Dirac-Frenkel a posteriori bound, familiar from variational
quantum time evolution); the claimed novelty is the *object* it is applied
to, a partition subspace of the computational basis, and what that buys
next.

### 4b. What leakage actually is: Theorem 3 (variance identity)

Theorem 3 is the conceptual heart. On a partition subspace the projection
residual has an exact statistical meaning: the per-class leakage equals a
**within-class variance**:

$$\lambda(v')^2 \;=\; \frac{1}{M_{v'}}\sum_{v} M_v\,\operatorname{Var}_{y\in
S_v}\!\big[g_{v'}(y)\big], \qquad g_{v'}(y)=\sum_d f_d(\beta)\,n(y;d,v').$$

Unpack it: $g_{v'}(y)$ is the **mixer-weighted neighborhood profile** of
bitstring $y$: "how strongly does the mixer connect $y$ back to cost-$v'$
states, summed over distances, weighted by $f_d(\beta)$." Leakage is large
exactly when **bitstrings that share a cost $v$ have *different* $g_{v'}$
values**, i.e. when "same cost" fails to imply "same neighborhood." That is
the precise, quantitative version of the proxy's founding intuition.

The crucial subtlety, and a genuinely useful result, is the **weighting**.
$\lvert f_d(\beta)\rvert = (\cos\beta)^{n-d}(\sin\beta)^d$ falls off fast
in $d$ for small $\beta$, so **only the small-distance (local) part of the
neighborhood variance matters** in the regime where the proxy is used.

> **Consequence for the G-RIPS triangle/Gaussian fitting:** fitting a shape
> to $N$ by *unweighted entrywise MSE* optimizes the **wrong norm.** It
> spends effort matching large-$d$ entries that $f_d(\beta)$ multiplies by
> nearly zero. This is the theoretical seed of the fitted-shape paradox and
> of the E017 direction experiment (section 5).

### 4c. Why it's so accurate at small angles: the Lemma and the cubic Corollary

This is the prettiest part and worth the five minutes.

**Step 1: half of every layer never leaks.** A QAOA layer is
$U = B(\beta)\Phi(\gamma)$. The phase separator $\Phi(\gamma)=e^{-i\gamma
C/2}$ is diagonal in cost, so it multiplies each class state $\ket{v}$ by a
single phase; it maps homogeneous states to homogeneous states **exactly,
on every graph, at every angle.** *All leakage comes from the mixer
$B(\beta)$.*

**Step 2: the mixer's first-order leak cancels, and that's a MaxCut
identity, not luck.** Expand the mixer to first order:
$B \approx I - i\beta\sum_j X_j$. The $\sum_j X_j$ part connects $y$ to its
$n$ single-bit-flip neighbors. To know how much this leaks, you need the
distribution of *neighbor costs* $c(y\oplus e_i)$. The **Lemma** computes
its first moment exactly:

$$\sum_{i=1}^n c(x\oplus e_i) = (n-4)\,c(x) + 2m \quad\text{for every } x.$$

The sum of all single-flip-neighbor costs is a **fixed linear function of
$c(x)$**, so it is *the same for every bitstring of a given cost*, i.e.
class-constant, i.e. **annihilated by $(I-P)$.** (Proof is one line:
flipping bit $i$ toggles the edges at vertex $i$; sum over $i$ using
$\sum_i\deg(i)=2m$. It's a consequence of MaxCut being **2-local**, not of
randomness. A two-step application of the same lemma also kills the
$O(\beta^2\gamma)$ double-flip channel, which the error form needs.)

**Step 3: the survivor is third order.** With the constant and linear
moments killed, the first term that leaks is the **second moment** of
neighbor costs. The result (**Corollary**):

$$\lambda_1 = \frac{|\beta|\,\gamma^2}{8}\sqrt{V_2}\,\big(1+O(\gamma)+O(\beta)\big),
\qquad V_2 = 2^{-n}\sum_v M_v\operatorname{Var}_{S_v}[s_2],\;
s_2(y)=\textstyle\sum_i c(y\oplus e_i)^2.$$

So compression error is **third order in the angles** ($\beta\gamma^2$),
with a prefactor that is a within-class variance of *squared* neighbor
costs. Measured exponents: $0.98$ in $\beta$, $1.92$ in $\gamma$ (E007).
This is the quantitative sharpening of Sud et al.'s own observation that
$k$-local QAOA states are Perfectly Homogeneous to leading order in
$\gamma$: now with the cancellation mechanism, the constant, and the
$\beta$-dependence.

**Bonus consequence, with its epistemic status flagged in the paper:** for
a linear-ramp schedule with fixed endpoints, each $\lambda_\ell$ depends
only on $(\gamma_\ell,\beta_\ell)\approx f(\ell/p)$, so
$\sum_\ell\lambda_\ell \approx p\cdot\overline{f}$: total leakage grows
**linearly in depth $p$.** The paper says plainly that this treats
$\lambda_\ell$ as a function of the layer angles alone (exact at $\ell=1$,
a small-leakage approximation afterwards), so the linear law is a
*prediction validated empirically*, not a corollary. Predicted ratio
$30/20 = 1.5$; measured $1.505$ on pooled sums, per-instance
$1.54\pm0.26$ SD over 70 index-paired runs (E011, E033).

### 4d. Opening the $V_2$ black box: the two codegree lemmas (E015)

The Corollary leaves $V_2$ as an opaque variance. Two exact,
machine-verified lemmas (to $10^{-10}$ on 280 instances; the cubic law
holds to $<0.3\%$ there) reduce it to pure graph structure:

- **Variance reduction.** Write $\delta_i(y) = c(y\oplus e_i) - c(y)$ for
  the single-flip cost changes and $T(y) = \sum_i \delta_i(y)^2$. Then
  $s_2 = (n-8)c^2 + 4mc + T$, and the first two terms are functions of the
  cost alone, so they vanish inside a within-class variance: **$V_2$ is
  just the within-class variance of $T$.** Only *squared single-flip cost
  changes* matter.
- **Codegree form.** In spin variables $\sigma_j=(-1)^{y_j}$, one line of
  algebra gives $\delta_i = \sigma_i\sum_{j\sim i}\sigma_j$, hence
  $T = \sum_{j,k} A_{jk}\sigma_j\sigma_k$ where $A_{jk}$ counts **common
  neighbors** of vertices $j,k$ (the codegree). Over uniform bitstrings
  this gives *exactly* $\mathbb{E}[T]=2m$ and
  $\mathrm{Var}(T) = 2\sum_{j\ne k}A_{jk}^2$, a closed form in $O(n^2+m)$,
  no $2^n$ enumeration. For $G(n,p_e)$ it evaluates to $\to 8p_e^2m^2$.

**The conditioning story (the mechanism behind the density law).** $V_2$
is the *cost-conditioned* version of $\mathrm{Var}(T)$. Unconditionally,
$\sqrt{\mathrm{Var}(T)}/m$ grows with density. But conditioning on the
cost removes a fraction that also grows with density; the two nearly
cancel, leaving the measured flat law $\sqrt{V_2}\propto m$. Empirically
the removed piece $\mathrm{Var}(T)-V_2$ tracks $\tau^2/m$ ($\tau$ =
triangle count) at Pearson $0.994$ (E015), but the triangle channel alone
undershoots badly; the next Proposition pins the correction. Slogan for the
paper: *"density drives compression error" sharpens to "squared codegrees,
conditioned on cost, drive it."*

### 4e. The quadratic conditioning bound (Proposition, E035)

This is the piece that turns the density law from an observation into a
formula. Everything runs through the **edge-spin sum** $S = m - 2c$, a
linear function of the cost, so conditioning on $c$ is the same as
conditioning on $S$. Four exact identities over uniform bitstrings (with
$\tau$ the triangle count and $c_4$ the 4-cycle count):

$$\mathbb{E}[S^3]=6\tau,\quad \mathrm{Var}(S^2)=2m^2-2m+24c_4,\quad
\mathrm{Cov}(T,S)=6\tau,\quad \mathrm{Cov}(T,S^2)=\mathrm{Var}(T).$$

The last one is the operative one: it is a **quadratic channel present even
in triangle-free graphs**, which is exactly why the crude linear estimate
($36\tau^2/m$) undershot the correction so badly (it captures only
18--64%, family-dependent).

**The Gram bound.** Project $T$ onto $\mathrm{span}\{S,\,S^2-m\}$. The
variance of that $L^2$ projection is $v^\top G^{-1} v$ with
$v = (6\tau,\ \mathrm{Var}(T))$ and $G$ the Gram matrix built from
$\mathrm{Var}(S)=m$, $\mathbb{E}[S^3]=6\tau$, and $\mathrm{Var}(S^2)$.
Since any projection onto a subspace is dominated by the projection onto
*all* functions of $S$, which is exactly
$\mathrm{Var}(\mathbb{E}[T\mid S]) = \mathrm{Var}(\mathbb{E}[T\mid c])$,

$$\mathrm{Var}\bigl(\mathbb{E}[T\mid c]\bigr) \ge v^\top G^{-1} v,
\qquad\text{equivalently}\qquad
V_2 \le \mathrm{Var}(T) - v^\top G^{-1} v.$$

**Proof idea.** Expand $T - 2m = 2\sum_{j<k}A_{jk}\sigma_j\sigma_k$ (the
codegree form) and $S = \sum_{(a,b)\in E}\sigma_a\sigma_b$; over uniform
spins **only even monomials survive**, so every moment is a graph count.
$\mathbb{E}[\sigma_j\sigma_k S]$ is an edge indicator, and summing
$A_{jk}$ over edges counts each **triangle at each of its three edges**,
giving $6\tau$ with the prefactor. $\mathbb{E}[\sigma_j\sigma_k S^2] =
2A_{jk}$ (**ordered two-edge paths** from $j$ to $k$), whose
$A_{jk}$-weighted sum is exactly $4\sum_{j<k}A_{jk}^2 = \mathrm{Var}(T)$.
$\mathbb{E}[S^3]$ counts ordered triangles; $\mathbb{E}[S^4]$ counts edge
pairings plus oriented 4-cycles.

**Why it matters practically:** every ingredient ($m$, $\tau$, $c_4$,
codegrees) is computable from the graph in **polynomial time**; no $2^n$
enumeration anywhere. Across 140 instances of the seven families the bound
captures 91--100% (mean 97%) of the exact conditioning correction, with
zero violations. For $G(n,p_e)$ it evaluates asymptotically to

$$V_2/\mathrm{Var}(T) \;\to\; 1 - \frac{p_e(1+4p_e-2p_e^2)}{1+6p_e^2-4p_e^3},$$

giving $0.375$ at $p_e=\tfrac12$ against $0.352$ measured at $n=12$--$14$.
The empirical density law becomes an explicit formula, up to a residual of
at most 10% (the degree-$\ge3$ residue: 3% on average, up to 9%). The
expansion also resolves the empirically observed "$\approx60\,\tau^2/m$" as
its finite-$n$ value; full identities, including the quartic split
$T = S^2 + m - R_4$, are in the companion note
`research/v2_conditioning_derivation.md`.

### 4f. Weighted MaxCut: the theory travels, the object doesn't (Remark, E037)

All of 4c--4e extends verbatim to integer-weighted MaxCut under mechanical
substitutions:

- **Neighbor-sum Lemma:** $m \to W$ (total weight).
- **Second moments:** replace $m$ by $\sum_e w_e^2$ and pick up a
  $\sum_e w_e^4$ term: $\mathbb{E}[T] = 2\sum_e w_e^2$,
  $\mathrm{Var}(S) = \sum_e w_e^2$, and
  $\mathrm{Var}(S^2) = 2\bigl(\sum_e w_e^2\bigr)^2 - 2\sum_e w_e^4 + 24c_4^w$.
- **Codegrees:** $A_{jk} \to \sum_i w_{ij}w_{ik}$; triangle and 4-cycle
  counts become the corresponding weight products.

All identities were machine re-verified to $10^{-9}$ on 140
integer-weighted instances (bound capture 92.9--99.6%, cubic law to 0.5%,
measured ratio $0.998\pm0.001$). The caveat that delimits the whole
enterprise: generic **continuous** weights collapse every cost class to a
bitstring-complement pair, so the compression itself becomes vacuous. Cost
degeneracy is a requirement of the *object*, not of the theory. (This is
precisely the degeneracy the quantile-binned classes of section 5 restore.)

### 4g. From state error to parameter error: the two-point regret certificate (E036, E038)

Everything so far bounds the **state** error, but parameter setting
consumes an **argmax**. The paper bridges the two with a certificate
evaluated at just two points (subsection "From state error to parameter
error", Proposition "Two-point regret certificate").

**Statement.** Let $F(\theta) = \langle\psi_p|C|\psi_p\rangle /
c_{\mathrm{opt}}$ be the true objective and $\hat F$ the proxy objective
computed from the **normalized** compressed state
$\chi = \phi_p/\lVert\phi_p\rVert$. If $\hat\theta$ maximizes $\hat F$ over
a candidate set containing a true maximizer $\theta^*$, then

$$\text{regret} \;=\; F(\theta^*) - F(\hat\theta) \;\le\;
\varepsilon(\theta^*) + \varepsilon(\hat\theta),$$

where for every $\theta$

$$\varepsilon(\theta) \le \frac{\lVert\psi-\chi\rVert}{c_{\mathrm{opt}}}
\Bigl(\sigma_\psi + \sqrt{\sigma_\chi^2 + \delta^2}\Bigr),
\qquad
\lVert\psi-\chi\rVert \le 2\sum_\ell \lambda_\ell(\theta),$$

with $\sigma_\psi, \sigma_\chi$ the cost standard deviations in the two
states and $\delta$ their mean-cost difference. Leakage is needed at only
the two relevant points, and both sit in the low-leakage corner in
practice.

**Proof idea (three moves).** (1) Decompose the regret into three terms,
$[F-\hat F](\theta^*) + [\hat F(\theta^*)-\hat F(\hat\theta)] +
[\hat F-F](\hat\theta)$; the **middle term is nonpositive** because
$\hat\theta$ maximizes $\hat F$. (2) For the two survivors, use the
identity (valid for normalized states and any scalar $a$)
$\langle\psi|C|\psi\rangle - \langle\chi|C|\chi\rangle =
\langle\psi-\chi|(C-a)\psi\rangle + \langle\chi|(C-a)(\psi-\chi)\rangle$;
**center at $a=\langle C\rangle_\psi$** and apply Cauchy-Schwarz, using
$\lVert(C-a)\chi\rVert^2 = \sigma_\chi^2+\delta^2$. (3) Bound the distance
via **Theorem 2**: $\lVert\psi-\phi\rVert \le \sum_\ell\lambda_\ell$, and
normalizing $\phi$ costs at most another
$1-\lVert\phi\rVert \le \lVert\psi-\phi\rVert$, hence the factor 2.

**The unnormalized variant.** The certificate as stated covers the
normalized proxy objective, which is also the convention recommended at
depth (section 5). For the unnormalized objective preferred at $p=1$ it
acquires one extra term,
$|a|(1-\lVert\phi\rVert^2)/c_{\mathrm{opt}} =
|a|\sum_\ell\lambda_\ell^2/c_{\mathrm{opt}}$ with
$a = \langle C\rangle_\psi$. (The paper marks this as an algebraic
addendum; since 2026-07-31 it is also machine-checked, by E036's
verify_addendum.jl at 714 points to 1e-10.)

**How to read it.** The certificate holds on all 280 instance-depth pairs
tested (machine-checked) but is honestly loose: median tightness is
$10$--$11\times$ over the actual regret, informative only where accumulated
leakage is a few percent. The slack has an exact anatomy (E038): with the
signed error field $e = F - \hat F$, regret decomposes exactly as
$[e(\theta^*) - e(\hat\theta)]$ minus the proxy's own margin
$[\hat F(\hat\theta) - \hat F(\theta^*)]$, and measurement shows three
stacked effects. The normalized proxy *overpredicts* essentially everywhere
($e<0$ at 557/560 individual argmax evaluations; both errors negative on
277/280 instance-depth pairs), so the two errors share a sign and partially
cancel (median cancellation factor $\approx0.6$); a **winner's curse**
concentrates the error at the proxy's own argmax
($|e(\hat\theta)| \approx 3.8$--$4.7\times |e(\theta^*)|$, because
maximizing $\hat F$ selects points of maximal overprediction); and the
margin absorbs roughly 40% of what remains. Regret is small because of this
stacked selection structure, which no pointwise norm, and no triangle
inequality, can see. Sharpening the certificate means modeling the
selection effect, the natural next theory question.

*(Note for the rewrite: the earlier draft's closing "why homogeneity
happens" ladder subsection is no longer a subsection of the QINP paper; its
content is distributed across the results above. If you want it as a
narrative device: the initial state is exactly homogeneous; the phase
separator preserves homogeneity exactly; the mixer is the only leak source
and only within-class variance of the weighted neighborhood enters; MaxCut
2-locality protects the first order on any graph; randomness buys the $V_2$
self-averaging rung, now pinned to within a few percent by E035; and good
schedules self-select the low-leakage corner.)*

---

## 5. What the experiments actually show (organized by question)

The paper's section 5 is a wall of numbers; here it is reorganized as the
questions each result answers. Setup, for context: seven graph families
(ER(0.5), ER(0.25), BA($k{=}2,4$), WS($k{=}4$; rewiring $0.1, 0.5$), random
3-regular), $n=12$--$20$, 20--30 instances per cell, fixed seeds. *Regret*
= grid ceiling minus proxy-set AR on the same schedule grid ($40\times40$
at $p=1$; $8^4$ linear-ramp endpoints at $p=3$). *Value-added* = proxy-set
AR minus the exact mean AR of random balanced partitions (already
0.75--0.76 on dense families; E002). Unless stated otherwise, proxy
argmaxes use the *unnormalized* objective (Sud et al.'s convention); the
paper recommends the normalized convention at depth. Family means carry
standard errors of 0.003--0.008 AR (value-added and AR means) and about
0.001 (regret cell means); only differences exceeding twice the relevant
error are interpreted.

**Q: Are the grid ceilings a fair reference? (E034)**
Yes, validated: Nelder-Mead refinement of the true objective from the grid
argmax raises the $p=1$ ceiling by at most 0.0007 AR (mean 0.0003) and the
$p=3$ ramp ceiling by 0.002 on average (max 0.005); releasing the
linear-ramp restriction entirely gains a further 0.007 on average (up to
0.05 on sparse instances). So $p=3$ regrets are against the best
coarse-grid linear ramp, within about 0.01 of the unrestricted continuous
optimum on average. The paper states this openly.

**Q: Is compression error about being Erdos-Renyi, or something simpler?**
**Density, full stop.** $\lambda/(\beta\gamma^2 m)$ is constant to
$\pm40\%$ across all seven families (E007). Surprisingly, ER(0.5), the
family the heuristic was built for, is the *worst*-compressing; random
3-regular the best. (Caveat the paper states: this fixes the angles;
practitioners often shrink $\gamma$ with density, which absorbs part of the
effect since $\lambda\propto\gamma^2 m$ is invariant under
$\gamma\mapsto\gamma/\sqrt{m}$.) Density is a known covariate of raw QAOA
*performance*; the claim here is the first density law for the error of
*compressing* the dynamics.

**Q: Is the telescoping bound actually usable, or just true?**
**Usable.** Slack 3--6$\times$ at $p=20$ over 840 runs, median about
$4\times$, never vacuous; near-functional predictor within each ramp
regime; linear-in-depth confirmed (1.505 vs 1.5). At fixed schedule,
accumulated leakage grows *sublinearly* in $m$ (+13--30% from $n=16$ to 20
while $m$ grows 30--58%), and small-ramp overlaps at $p=30$ stay
0.71--0.81 across $n=16$--$20$ (E003, E011).

**Q: When the proxy degrades, is it because the state got complicated, or
because the fixed cost-class frame points the wrong way?**
**Wrong way (mis-aim), not complexity.** Stacking $p=20$ trajectory states
and doing PCA shows effective dimension 2--4 (99% of energy), while the
$(m{+}1)$-dim cost-class subspace captures 79--99% depending on angles; a
2--3-dimensional trajectory-PCA subspace matches the entire cost-class
frame (E009). The state stays simple; the *fixed* frame drifts off the tiny
moving subspace. (Trajectory PCA is a diagnostic benchmark, not a method;
it needs the states the proxy exists to avoid computing.)

**Q: Does low compression error mean good parameters?**
**At $p=1$, yes, and the ranking holds with pre-committed statistics.**
Exact-compression regret is 0.028--0.038 at $n=12$--$14$ and 0.032--0.055
at $n=16$--$18$, roughly doubling at $p=3$ ramps (0.06--0.11), with
positive family-level value-added throughout (0.02--0.13, dense families
lowest; under sampled-$N$ one $G(16,0.5)$ instance dips to a marginal
$-0.007$). Regret grows slowly and *monotonically* with $n$: +0.005 (dense
ER) to +0.018 (3-regular) from $n=12$ to 18 at $p=1$, small but resolvable
against the $\sim0.001$ per-cell SEs. (The paper explicitly corrects its
own earlier "flat in $n$" reading.) Families ranked by leakage rank by
regret at $p=1$: Spearman $\rho=0.96$ at $n=12$ and $0.86$ at $n=14$, with
instance-bootstrap 95% CIs $[0.64,0.96]$ and $[0.61,0.93]$;
$\rho\gtrsim0.8$ was a pre-committed acceptance criterion, not post-hoc
(E002, E004, E010, E033). **Honest limit:** at deeper $p$ and larger $n$
the family regret differences compress to within a few 0.01 and the
ranking signal fades with them; a scope limit, not a sign reversal.

**Q: At depth, is regret a fidelity problem at all? (E014, E033, E040)**
**No. It is an argmax-transfer problem, and fidelity decouples.** Regret
tracks the proxy's **argmax displacement** from the true optimum at both
depths: pooled Spearman $\rho=0.74$ (95% CI 0.64--0.82) at $p=1$ and 0.65
(0.54--0.74) at $p=3$, essentially unchanged under family-and-size
demeaning (0.63, 0.67), within-cell averaging, and cluster bootstrap. The
state-fidelity deficit is **family-confounded both ways**: pooled
$\rho=0.39$ at $p=1$ drops to 0.07 within family-size cells, while the
$p=3$ pooled value of $-0.02$ hides a moderate within-cell 0.46
(proxy-chosen schedules compress fidelity into a narrow band across
families). Landscape flat-peak robustness predicts nothing
($|\rho|\le0.24$). At $n=16$ (35 instances) the displacement-regret
correlations are 0.84 and 0.76 under the normalized-objective regret,
though only 0.25 and 0.47 under the unnormalized convention; the paper
states the geometry is stable one size up *for the depth-recommended
convention*, not uniformly. Punchline: the leakage calculus bounds the
*state* error while parameter quality lives in *parameter space*, which is
exactly why a lower-fidelity model with a better-placed argmax (the
analytical $N$ off-ER) beats the exact compression.

**Q (model error): how bad is using the analytical $N$ off-ER?**
**Its argmax is excellent, its values are garbage.** The one-line thesis of
the paper's model-error subsection: *across every model class tested, only
the argmax transfers; no value-, norm-, or fit-based quantity certifies
it.* The analytical formula (with an effective edge probability
$2m/n(n{-}1)$) gets regret 0.005--0.02 at $p=1$ on BA/WS/3-regular/sparse
ER, *better* than exact compression, because the smooth model regularizes
instance noise. But on **dense ER(0.5)** its argmax lands on an unphysical
spurious peak ($\langle C\rangle=93$ on a 38-edge graph, squared state norm
inflated $7.5\times$). The norm certificate *detects* this but **cannot
repair it**: on sparse families the model's calibration breaks everywhere
while its argmax stays excellent, so rejecting norm-inflated or physically
impossible predictions rejects the answer along with the artifact (regret
rising from about 0.01 to 0.03--0.31 across families and depths). Reported
as a negative result (E004, E005, E006). One refinement from E017:
"absolute" is load-bearing; the raw value *landscape* is
argmax-informative, and post-hoc renormalization discards that signal
(dividing the analytical model's predictions by its state weight triples
its regret, $0.05\to0.16$).

**Q: does fitting better shapes to $N$ help? (E012 + E017, the paradox
resolved)**
**No, and now we know exactly why.** Fitting a triangle/Gaussian to $N$ by
entrywise MSE with the original objectives is harmful (triangle: mean
$p{=}1$ regret $0.054\to0.211$, worse on 136/140) or inert (Gaussian:
tenfold MSE improvement, argmax unchanged on 140/140), and across the model
zoo **no scalar mismatch norm predicts regret**; only argmax displacement
does ($\rho\approx0.7$). E017's controlled perturbations supply the
mechanism: corrupting the exact $N$ at *matched* entrywise magnitude but
controlled direction, the $f_d(\beta)$-weighted transfer error of Theorem 3
governs regret wherever it has dynamic range (pooled $\rho=0.74$ vs 0.38
for MSE). Directions invisible to the weighted norm are free even at 50%
entrywise error; aligned directions at 1% already move the argmax.
Entrywise MSE doesn't just fail to predict, it **misorders** the two shape
families against regret in 105/150 instances (worse than chance, the
G-RIPS paradox exactly). But the weighted norm also explains E012's null:
every realistic shape fit *saturates* it (relative weighted error
$\ge0.92$ even when fitted directly in that norm), so refitting in the
correct norm does not rescue the shapes. **The shape-fitting program fails
on representational grounds, not from a wrong objective**, which is one
more reason the sampled-$N$ recipe needs no shape class. Two byproducts:
Theorem 1 contractivity held on all 150 instances (exact-$N$ state weight
$\le0.984$), so any norm inflation in the wild certifies model error; and
normalization is an *instrument, not a recipe*.

**Q: what is the best proxy variant, where any is competitive?**
**Sampled $N$ ($S\approx10$ stratified samples per cost class,
$O(S\,m\,2^n)$ instead of $O(4^n)$) plus the empirical cost distribution.**
It matches exact-$N$ parameter setting within 0.011 AR at $p=1$ and 0.016
at $p=3$ (pooled; SEs $\le0.004$), with pooled regret *lower* by
0.003--0.006 (lower in 26 of 28 family-size cells): sampling regularizes
like the analytical model, without its dense-graph artifact, on every
family (E010, E033). Leakage itself is estimable the same way: $S=5$
reproduces the aggregate leakage fraction to median 3.2% error (E008).
Honest scope, stated in the paper: the cost is still exponential, so it
serves where statevector passes are affordable but $O(4^n)$ is not (roughly
$n\lesssim30$); whether $N$ admits a polynomial-time Monte Carlo estimator
is open.

**Q: why normalize the proxy objective at depth? (E036, E039, E040)**
**Because the unnormalized objective carries a leakage bias.** Along the
depth grid the compressed norm varies substantially from schedule to
schedule, so the unnormalized objective conflates "high predicted value"
with "low leakage" and favors small-angle schedules. Dividing by the
tracked compressed norm (free, from Theorem 2's bookkeeping) removes the
bias exactly where leakage accumulates: pooled $p=3$ regret drops
$0.080\to0.044$ (better on 134/140); at $p=1$ the unnormalized convention
keeps a mild edge (0.031 vs 0.047). Replicates at $n=16$:
$0.090\to0.060$, better on 32/35 (E040). Rule: **normalize at depth, keep
unnormalized at $p=1$.** The paper connects this to the
instrument-not-recipe caution: normalization divides by a norm the exact
compression tracks *correctly*; dividing an analytical model's predictions
by its miscalibrated state weight discards argmax signal. Pushing further
fails (E039, negative): the residual overprediction is only weakly
correlated with norm loss (per-instance $\rho$ from $-0.3$ to 0.6), and a
pooled linear norm-loss correction worsens regret on every instance it
moves. Normalization removes the only norm-visible part of the bias; the
winner's-curse remainder is a selection effect no pointwise correction can
see. (Third independent confirmation, after the fitted shapes and the
filters, that calibrating the proxy's *values* harms its *argmax*.)

**Q: how does the proxy compare to plain parameter transfer? (E018, E034)**
**Transfer is at least as good, twice replicated with different source
constructions.** (1) E018: angles grid-optimized on a single $n=10$
same-family source achieve pooled regret 0.007 (median 0.003), modestly
below every proxy variant (0.03--0.05; both methods near the ceiling, so
the gap is a few 0.01 of AR). The paper reads this as *agreeing* with Sud
et al.'s statistical tie at $p\le3$: with a single cheap source rather than
their robust median, transfer lands marginally closer, but neither
dominates. (2) E034: averaging the true grid-argmax angles of ten
brute-forced ER(0.5) $n=12$ instances and applying them verbatim to every
instance, family, and size gives pooled regret about 0.014 at $p=1$ and
0.006 at $p=3$, better than the exact compression on 134/140 instances at
$p=1$ and on all 140 at $p=3$. This is parameter concentration made
concrete: within a family the argmax barely moves, so any cheap view of a
source landscape suffices, and the proxy's per-instance distribution buys
no edge.

**Q: what about the two regimes those tests didn't cover, high depth and
beyond the exact-$N$ wall? (E019, the harsh result)**
**Both go to transfer, with no asterisks.** The paper runs the high-depth
head-to-head Sud et al. could not: their "transfer fails at $p=20$"
reflected that no 4-parameter ramp-transfer schedule existed in 2022, not
an inferiority of ramp transfer, which the intervening linear-ramp
literature established. With it, at $p=10, 20$ ($n=14$--$16$): their
monotonically increasing *absolute* AR stays intact, but measured as regret
against the best $p=20$ ramp, **every proxy variant collapses** (exact
compression 0.07--0.14, sampled-$N$ 0.11--0.21, analytical 0.18--0.37,
i.e. regret 0.07--0.37 across variants) while ramp transfer holds
$\le0.03$ (best-of-3 sources $\le0.007$ in every cell). Even the
zero-model-error compression misplaces the argmax by about 0.1 AR at
fidelity 0.7--0.8: the argmax-transfer failure mode in pure form. Beyond
the exact-$N$ wall ($n=22, 26$), transfer still wins every cell. One
genuine positive for the analytical model: on sparse regular graphs its
argmax *improves* with $n$ (regret $0.019\to0.017$ from $n=22$ to 26; its
class-level derivation is asymptotic), and its native binomial cost
distribution is load-bearing: substituting the empirical one degrades
3-regular $p=3$ regret from 0.05 to 0.27, its conventions mattering just as
with the fitted shapes.

**Q: does the calculus survive continuous weights? (E020, the binned proxy)**
**Yes: Theorems 2--3 hold for any fixed partition.** With i.i.d. $U[0,1]$
weights every cost class is a singleton, but taking $K$ equal-population
(quantile) bins of the weighted cost as the classes, the measured one-layer
leakage decomposes as
$\lambda^2\approx\lambda_{\mathrm{struct}}^2+O(1/K)$, a structural
compression floor plus a binning variance, and the $K$-dimensional binned
proxy (bin-label distributions through the same machinery, bin-mean
phases) recovers integer-cost regret levels by $K\approx64$--$128$: 0.038
vs about 0.032 unweighted on $G(16,0.5)$ at $p=1$; 0.103 vs 0.105 on
3-regular at $p=3$. Cost: $O(K^2n)$ per layer per schedule, independent of
$m$, so the weighted extension is *cheaper* to run than the integer proxy
on dense graphs.

**Q: do weights finally give the proxy a winning regime? (E021, E022--E026,
and the E032 audit; the retraction story)**
**No, and the paper corrects itself in the open.** Mild weights first:
with i.i.d. uniform or exponential weights, single-source transfer regret
rises three- to tenfold yet still beats the binned proxies in every cell
(E021, E022); only the zero-knowledge pooled schedule collapses on dense ER
at $p=3$. The July line then built an apparent positive result (E022--E026):
under heavy-tailed weights (Pareto $\alpha=1.5$) the binned proxy seemed to
win $p=1$ cells by 2--7$\times$, with a pre-registered dispersion boundary
at CV$^2\approx1$ (E023 located it at the infinite-variance point
$\alpha\approx2$; E025 showed via lognormal weights that finite-variance
dispersion suffices; E026 confirmed a pre-registered flip at
$\alpha\approx2.4$ on dense ER; E024 showed it surviving to $n=22, 26$).
**E032 audited the baseline and the niche dissolved.** The transfer
baseline had rescaled its $\gamma$-grid by each instance's *mean edge
weight*, which is tail-dominated under heavy tails, exactly the regime
tested. Rescaled instead by the cost-distribution width (the natural phase
scale, and what the scale-free quantile bins implicitly use), transfer
*ties* the proxy single-source and *wins* best-of-three on the same
$G(16,\tfrac12)$ Pareto instances: proxy 0.012; mean-rescaled transfer
0.090; cost-width-rescaled 0.015 single / 0.002 best-of-three. The
CV$^2\approx1$ "boundary" is near-tautological: it marks where the mean
stops being a usable scale, not an intrinsic QAOA crossover. The paper
collapses the whole arc to one honest paragraph: the binned proxy is
*scale-free* and therefore **matches well-scaled transfer without needing a
scale estimator: a convenience, not an advantage.** "An earlier version of
this experiment mistook the mean-rescaling handicap for a genuine
crossover; we correct it in the open."

**Q: is the deep high-dispersion corner intrinsically hard? (E028, E029)**
No: it's a **grid artifact**. E028 found every method bad (regret up to
about 0.3) at (high dispersion $\times$ $p=10$--$20$ ramps). E029 shows
compass refinement of the true objective from *any* grid start (proxy,
transfer, universal, or the grid ceiling itself) converges to essentially
one optimum: continuous ramp endpoints cut regret to 0.010--0.037 and
freeing all $2p$ angles finishes the job (residual $<0.006$ from every
start). The landscape is not trap-riddled, it is *sharp*: optima fall
between the endpoint-grid points, and a few thousand true evaluations erase
the corner.

**Q: can the proxy at least refine for free where evaluations are precious?
(E030, E031)**
**Yes, within the schedule family and with a guard.** Polishing the four
ramp endpoints on the (binned sampled) proxy's own landscape lifts true AR
by +0.05--0.12 for one confirmation evaluation (lognormal
$0.799\to0.915$; Pareto $0.837\to0.884$). Freeing all $2p$ angles on the
surrogate *overfits its model error* under extreme tails (Pareto full-$2p$
falls to 0.829, below its own ramp result). The composed recipe (map-init,
surrogate ramp-polish, two true evaluations, keep the better of initializer
and polished) never hurts and gains exactly where the surrogate has signal;
without the guard, polish *degrades* unweighted inits already at the
ceiling ($0.897\to0.837$). Paper statement: "the surrogate refines for
free ... if one stays within the schedule family and keeps the better of
initializer and polished."

**Q: is any of this MaxCut-specific? (E027)**
**Only the small-angle cancellation, as the theory says.** A Max-3-XOR
replication (Sud et al.'s second family; integer costs, native classes)
confirms the map: exact-compression $p=1$ regret 0.031--0.040 (MaxCut
magnitude); sampled-$N$ again matches or beats exact (7 of 8 cells);
transfer again dominant on unweighted instances (near 0.000 regret in most
cells; concentration is even stronger than MaxCut's); small-angle leakage
doubles when clause density doubles ($0.045\to0.092$ at $n=14$: the
density law again). Large-angle leakage is much higher than MaxCut's
(0.84--0.95 at $\gamma=1.0$), consistent with 3-local phases lacking the
MaxCut cancellation, the one genuinely MaxCut-specific piece.

**Q: what does the pipeline actually cost? (E013)**
The honest accounting that bounds the pitch. On one $G(20,0.5)$ instance
($m=103$; A100 + 8 CPU threads): exact $N$ 21.8 s ($O(4^n)$, about
$16\times$ per two qubits, hours by $n=24$); sampled $N$ ($S=10$) 0.44 s
($50\times$ cheaper); proxy sweep $p=1$, 1600 schedules, 0.38 s; $p=20$,
4096 schedules, 13.4 s; one statevector $\langle C\rangle$ at $p=1$ on GPU
0.3 ms, so the brute-force grid ceiling (1600 evaluations) costs 0.48 s,
no more than the proxy sweep *before* counting the cost of $N$. At
simulable sizes the proxy buys no wall-clock advantage; its cost regimes
are beyond-simulation sizes with the analytical $N$, and GPU-less
environments. Hardware shots do not separate methods: every offline
surrogate, transfer included, avoids them equally.

**Q: so what should a practitioner actually do?**
The paper's practical section now composes into this map (it is threaded
through 5.4/5.5 rather than stated as a numbered rule):

1. **Where a same-family solvable source exists** (every random family
   tested, at simulable sizes): **transfer** its angles; with weights,
   rescale the schedule by a robust cost scale (cost-distribution width,
   not mean weight). Transfer matched or beat every proxy variant tested,
   at $p=1$ through $p=20$, unweighted and weighted alike (E018, E019,
   E021, E032, E034).
2. **Among proxy variants, at low depth ($p\le3$), where statevector passes
   are affordable but $O(4^n)$ is not:** sampled $N$ ($S\approx10$) + the
   empirical cost distribution, normalized objective at depth. Matches
   exact within 0.011--0.016 AR with slightly lower regret; family-
   agnostic; no dense-graph artifact (E010). For continuous weights, the
   quantile-binned version with $K\approx100$ bins at $O(K^2n)$ (E020). At
   $p=10$--$20$ every proxy variant collapses (E019); don't use one there.
3. **Free refinement:** polish ramp endpoints on the surrogate's landscape,
   confirm with one true evaluation, keep the better of init and polished
   (E030, E031).
4. **Beyond classical reach:** only the analytical $N$ exists. Trust its
   argmax, never its values, and distrust it on dense Erdos-Renyi
   (E004--E006, E019).

**Q: stepping back, when does the proxy work at all?**
The paper's closing synthesis (Conclusion): the compression is faithful
wherever density-controlled leakage is small at the schedules worth
running; parameter quality is an argmax-transfer matter that decouples from
state fidelity at depth; model error lives in the error's *direction*
against the mixer-weighted norm. Calibrated against the strongest cheap
baselines, the proxy loses every homogeneous regime and does not surpass
transfer in any regime with a fair baseline; this *sharpens, not
overturns*, Sud et al.'s own report of parity. Its practical residue is
modest and honest: a scale-free per-instance surrogate that matches (not
beats) well-tuned transfer, free guarded polishing, the analytical model
beyond simulability, and, most durably, the error calculus itself.

---

## 6. What is proved vs. measured (so you can defend it)

- **Proved (exact, machine-verified):** Theorem 1 (proxy = compression, to
  $2.5\times10^{-16}$), Theorem 2 (telescoping bound + norm identity),
  Theorem 3 (variance identity), the neighbor-sum Lemma, the Corollary's
  cancellation (first-order leak identically zero, including the
  $O(\beta^2\gamma)$ double-flip channel), the two $V_2$ lemmas, the four
  moment identities and the quadratic conditioning bound (E035), their
  integer-weighted extensions (E037), the two-point regret certificate
  (E036), and the exact signed regret decomposition (E038's identity).
  These hold for *any* graph; no randomness needed. Theorems 2--3 hold for
  *any* partition of *any* cost function (which is what licenses the binned
  weighted proxy and the Max-3-XOR replication); only the cubic
  cancellation is MaxCut's own.
- **Measured, not proved:** the density law's self-averaging (now pinned to
  within a few percent by E035, but the residual and a rigorous
  self-averaging statement remain); linear-in-depth accumulation (a
  prediction validated, not a corollary); all regret numbers; family
  rankings; the sampled-$N$ recipe's quality; every transfer comparison.
- **The honest seams** (the paper states all of these):
  - The leakage-to-regret ranking is a 7-point Spearman resolvable only at
    $p=1$/small $n$; it fades at scale.
  - Experiments are unweighted MaxCut at $n\le20$ for the leakage anatomy;
    parameter-quality comparisons reach $n=26$ and $p=20$. The weighted
    results span uniform, exponential, and Pareto weights at $n=16$, 10
    instances per cell, and rest on one robust choice of transfer scale;
    the binning term's constant is empirical. Refinement and recipe checks
    use 5 instances per cell.
  - All comparisons are noiseless (offline angle selection); schedule
    robustness under hardware noise is not studied.
  - Trajectory PCA is a diagnostic, not a method.
  - The sampled-$N$ recipe's cost is still exponential.
  - The certificate is honestly loose (median 10--11$\times$), and the
    $n=16$ stability of the displacement geometry holds under the
    normalized convention only.

---

## 7. The open doors

1. **Instance-adapted low-rank frames: the cheap version fails (E016, in
   the Discussion).** The trajectory really is about 4-dimensional (an
   oracle rank-4 PCA frame captures about 0.99 of a depth-20 trajectory),
   so a good instance-adapted frame would beat the $(m{+}1)$-dim cost-class
   frame handily. But a rank-4 frame built from the first five layers is
   nearly orthogonal to the late-time subspace (largest principal angle
   71--88 degrees, growing with ramp) and captures *less* than the
   zero-cost cost-class frame. The low-rank subspace **rotates** across
   depth; any frame that beats the compression must model that rotation
   explicitly. A measured obstacle, reported as future-work guidance.
2. **The $V_2$ residual.** The quadratic bound captures 91--100% (mean
   97%); what remains open is the degree-$\ge3$ residue (3% on average, up
   to 9%) and a rigorous self-averaging statement, which would turn the
   density law fully into a theorem.
3. **A polynomial-time estimator of $N$** (or of leakage) would extend the
   sampled-$N$ recipe beyond the classically simulable regime; nothing
   rules it out.
4. **Sharpening the certificate.** The 10--11$\times$ looseness IS the
   argmax-transfer story quantified: overprediction with a shared sign,
   winner's curse at the proxy's argmax, margin absorption. A sharp regret
   theory would model the selection effect, not the pointwise error.
5. **Structured-instance transfer.** Every transfer victory is on
   concentrating random ensembles; whether genuinely heterogeneous or
   structured instances (outside any family with concentrating angles) give
   the per-instance surrogate a real niche remains untested at depth (the
   heavy-tail candidate niche died in the E032 audit; the deep-dispersion
   corner proved a grid artifact in E029).

---

## 8. The noise-regularization arc (E041--E043; folded into the paper 2026-07-31)

*Status: E041's README quarantined this claim until the mechanism was
resolved; E042/E043 resolved it, and the arc entered qinp_paper.tex on
2026-07-31 as one scoped paragraph at the end of Sec. 6.1 (best proxy
variant at $p\le3$; within a family, transfer still wins).*

- **E041 (sampled + normalized composes, and overshoots).** The two
  separately-validated ingredients combine better than either alone:
  sampled $N$ ($S=10$) + normalized objective is the best of all four
  {exact, sampled} $\times$ {raw, normalized} combinations at BOTH depths,
  beating even exact $N$ + normalized: $p=1$ mean regret 0.021 vs 0.047
  (123/7 win/loss, sign $p=8\times10^{-29}$); $p=3$ 0.024 vs 0.044
  (112/23), uniformly across all 7 families. With sampled $N$ the
  normalized objective helps at $p=1$ too, unlike with exact $N$. A proxy
  driven by *less data* sets *better* parameters than the exact
  compression.
- **E042 (mechanism: noise, not estimator bias).** Replicate-averaging the
  sampled $N$ (noise suppressed, bias kept) loses the gain (regret back
  near exact-$N$ level: 0.036/0.037 vs exact 0.042/0.040); the $S$-sweep
  is monotone with *smaller* $S$ better ($S=3$: 0.012/0.016; converging to
  exact by $S=300$); and noise does nothing for the raw objective at any
  $S$. Sampling noise disrupts the winner's-curse selection effect (E038)
  that the normalized proxy's systematic errors feed.
- **E043 (the noise must be structured).** Depth dissociates generic from
  structured noise: at $p=1$, i.i.d. multiplicative noise on exact $N$
  reproduces and even beats the sampled gain ($\sigma=0.05$--$0.1$: regret
  0.011--0.019, U-shaped in $\sigma$), but at $p=3$ i.i.d. noise only
  hurts (0.036--0.126 at every $\sigma$, vs sampled 0.023). The sampler's
  rows are genuine member profiles with exact marginals (each row sums to
  $\binom{n}{d}$), so its noise lives inside the physical constraint set;
  i.i.d. noise leaves it, and destroys the raw objective everywhere
  (0.14--0.30). At depth the gain requires the sampler's *structured*
  noise.

If folded into the paper, this would upgrade the practitioner recipe
(sampled + normalized at both depths, with smaller $S$ apparently better)
and add a genuinely novel finding: **structured subsampling noise as a
regularizer against the proxy's winner's curse.** It would also interact
with the paper's current "normalize at depth, raw at $p=1$" rule, which is
an exact-$N$ statement. Caveats before promoting: E042/E043 are $n=12$
only; one noise model; single sampling budget landscape mostly unexplored.

---

## 9. Symbol cheat-sheet (friendlier than the paper's Table 1)

| Symbol | Say it as | Meaning |
|---|---|---|
| $n,\,m,\,p$ | "qubits, edges, depth" | vertices; edges (= max cost); QAOA layers |
| $c(x)$ | "cost of $x$" | cut value of bitstring $x$ |
| $S_v,\,M_v$ | "cost class $v$, its size" | all $x$ with $c(x)=v$; how many there are |
| $\ket{v}$ | "class state" | uniform superposition over $S_v$ (normalized) |
| $\mathcal{H}_{\mathrm{hom}}$ | "homogeneous subspace" | states constant on each cost class ($\dim\le m{+}1$) |
| $P$ | "the blur" / "class averaging" | projector onto $\mathcal{H}_{\mathrm{hom}}$: average within each bucket |
| $U_\ell=B(\beta)\Phi(\gamma)$ | "a QAOA layer" | mixer times phase separator |
| $\Phi(\gamma)$ | "phase separator" | diagonal; **preserves homogeneity exactly** |
| $B(\beta)$ | "mixer" | the **only** source of leakage |
| $n(x;d,v)$ | "neighborhood profile" | # bitstrings at distance $d$ from $x$ with cost $v$ |
| $N(v';d,v)$ | "homogeneous distribution" | class-average of $n$; the recursion's coefficients |
| $Q_\ell(v)$ | "proxy amplitude" | the one number per cost the proxy tracks = blurred amplitude |
| $\lambda_\ell$ | "leakage" | norm that leaves $\mathcal{H}_{\mathrm{hom}}$ in layer $\ell$ |
| $\phi_\ell$ vs $\psi_\ell$ | "compressed vs true" | proxy (= blurred) trajectory vs real QAOA |
| regret | "regret" | (best AR on the grid) $-$ (AR at the proxy's chosen angles) |
| value-added | "value-added" | (proxy AR) $-$ (random balanced-partition AR $\approx0.75$) |
| argmax transfer | "does the peak land right?" | the only thing parameter setting actually consumes |
| $\delta_i(y)$ | "single-flip change" | $c(y\oplus e_i)-c(y)$; spin form $\sigma_i\sum_{j\sim i}\sigma_j$ |
| $T(y)$ | "flip-energy sum" | $\sum_i \delta_i^2$; $V_2$ = its within-class variance |
| $A_{jk}$ | "codegree" | # common neighbors of $j,k$; $\mathrm{Var}(T)=2\sum_{j\ne k}A_{jk}^2$ |
| $S = m-2c$ | "edge-spin sum" | conditioning on $c$ = conditioning on $S$; carrier of the E035 bound |
| $\tau,\,c_4$ | "triangles, 4-cycles" | the graph counts in the conditioning bound |
| $e(\theta)=F-\hat F$ | "signed error field" | overprediction is $e<0$; regret $= \Delta e -$ margin exactly |
| $K$ | "bins" | quantile cost classes for continuous weights; proxy cost $O(K^2n)$ |

---

## 10. The paper's skeleton (for rewriting it from scratch)

The canonical draft (`qinp_paper.tex`, single-column, ~20 pp) is organized
so that each section owns one move. If you rewrite, this is the
load-bearing order:

**The abstract** is shaped thesis-first, no suspense: (i) the proxy is an
*exact orthogonal compression*, not a substitution heuristic; (ii) the
error calculus it yields (leakage at $O(2^n)$, the telescoping bound with
measured slack 3--6$\times$, the variance identity naming the norm in which
distribution models must be judged); (iii) the MaxCut small-angle law
($O(\beta\gamma^2)$, linear in depth, density-governed); (iv) parameter
quality is argmax transfer, decoupling from fidelity by $p=3$, which is why
MSE fitting backfires and no scalar mismatch certifies a model; (v) the
honest calibration: the proxy is *not* a faster parameter setter, plain
within-family transfer is at least as good as every variant, extending Sud
et al.'s low-depth parity to high depth; the lasting contribution is the
leakage certificate.

1. **Introduction.** The hook is the two open questions (what is the proxy
   / when does it work), quoted against Sud et al.'s own "no longer
   unitary, analogues of amplitudes" language. Explicitly a
   characterization-with-error-calculus paper, not a speedup paper (the
   GPU evaluates the true objective about as fast as the proxy scores it;
   all offline methods are equally shot-free). Contributions list:
   exactness (verified to $2.5\times10^{-16}$); error calculus (slack
   3--6$\times$, never vacuous); small-angle structure (exponents 0.98 /
   1.92; linear depth 1.505 vs 1.5); experimental anatomy (density not
   family; regret 0.03--0.05 at $p=1$ growing mildly in $n$; effective
   dimension 2--4); practical findings positive and negative (analytical
   argmax robust / values decalibrated; filters fail; sampled $N$ best
   proxy variant at $p\le3$; every variant collapses at $p=10$--$20$;
   transfer matches or beats throughout, weighted included via robust
   rescaling; quantile binning with $\lambda^2\approx
   \lambda_{\mathrm{struct}}^2+O(1/K)$ and ~100 bins). Footnote: all
   claims trace to `research/experiments/`, cited as E001--E043 with gaps.
2. **Background and related work.** Short QAOA and proxy subsections; the
   proxy subsection now states Sud et al.'s two claims precisely
   (statistical tie with median transfer at $p\le3$; monotone ramp AR to
   $p=20$ with no transfer table then existing), which the paper "takes as
   its starting point." Related work does the novelty positioning: exact
   symmetry reductions and quotient-graph walks as the exact special case
   (Grover mixer = zero leakage, transverse field makes $N$ the relevant
   object); lumpability / model order reduction as the classical
   mathematics; constrained bisimulations as exact reductions vs our
   approximate error calculus; Dirac-Frenkel a posteriori bounds as the
   proof lineage (novelty claimed only for the object); parameter-setting
   toolbox (transfer, ramps, concentration) with the note that grid
   ceilings, not transfer, are the regret reference, so regrets are
   conservative; classical surrogates of QAOA as the proxy's family;
   Krueger-Mauerer as the closest landscape-surrogate neighbor (our
   contribution on that axis is diagnostic: *which* error norms control
   argmax quality). A glossary table (Table 1) keys the vocabulary.
3. **Section 3: the proxy is an exact subspace compression.** Setup
   ($S_v$, $M_v$, class states, $P$, empirical $N$); Theorem 1 with the
   two-line induction proof and transfer-matrix form; three delimiting
   remarks (same-instance $N$ only, model error $E$ named here; unattained
   costs inert; do not renormalize); Proposition 1 (norm inflation is a
   zero-cost model-error certificate) with the pointer that its converse
   fails instructively.
4. **Section 4: the error calculus of leakage.** Theorem 2 (telescoping +
   Pythagoras) and the instrument paragraph with the Dirac-Frenkel
   attribution; Theorem 3 (variance identity + the wrong-norm warning);
   subsection "Small-angle structure for MaxCut": neighbor-sum Lemma,
   cubic Corollary (with the epistemic flag on the linear-depth step), the
   two $V_2$ lemmas, the conditioning Remark, the quadratic conditioning
   bound (Proposition, E035) with the ER asymptotic and the companion-note
   footnote, the weighted-MaxCut Remark (E037); subsection "From state
   error to parameter error": the two-point regret certificate
   (Proposition, E036), its measured tightness and the E038 signed
   anatomy in prose, and the unnormalized-variant addendum. (No "why
   homogeneity happens" ladder subsection in this version.)
5. **Section 5: experimental anatomy.** 5.1 Setup: seven families, sizes,
   grids, regret/value-added definitions, ceiling validation (E034),
   objective convention (unnormalized unless stated; E036), verification
   levels, SE discipline. 5.2 Compression error: leakage maps + density
   law (Fig: anatomy), bound tightness + depth scaling (Fig: bound),
   sublinear-in-$m$, trajectory PCA. 5.3 Compression error and parameter
   quality: regret table (Table 2), monotone-$n$ growth, leakage-regret
   ranking with CIs and scope (Fig: ranking), argmax displacement vs
   fidelity with all controls and the $n=16$ stability sentence (Fig:
   argmax). 5.4 Model error: thesis line ("only the argmax transfers"),
   analytical robustness + dense-ER artifact, filter negatives, fitted-
   shape paradox dissected by E012 + E017 (Fig: directions), the
   instrument-not-recipe byproducts. 5.5 "What should a practitioner do?":
   sampled-$N$ recipe with honest scope; the normalization rule (E036,
   E039, E040); the transfer calibration (E018 + E034) and Sud
   reconciliation; the high-depth head-to-head (E019); the weighted
   paragraph: binned proxy (E020, Fig: binned), mild-weight transfer
   (E021), the scale-audit correction told openly (E018/E022/E032);
   refinement findings (E029; E030/E031 guarded polish). 5.6 "What the
   pipeline costs": timing table (Table 3) and the bounded pitch.
6. **Discussion.** "What the compression view buys" (conserved quantities,
   cheap error, provable structure, the compression/model split shown
   non-academic); "Limits" (scopes and instance counts per result family;
   Max-3-XOR replication E027 scoping what is MaxCut-specific; noiseless
   comparisons; the E016 rotating-subspace obstacle; the remaining
   $V_2$ residue).
7. **Conclusion.** The candid map, ending on the durable contribution: the
   error calculus, "a way to know, before trusting any cost-class
   surrogate, how much of the state it is silently discarding."

Style rules the draft follows: every quantitative sentence carries a
`% E0xx` comment; negative results and self-corrections are reported as
findings, in the open; each theorem's scope (any graph / any partition vs
MaxCut-specific vs random-like) is stated where it is used.

---

## 11. Key numbers, with provenance (post-reconciliation, 2026-07-31)

| Claim | Number | Experiment |
|---|---|---|
| Thm 1 exactness verified | $2.5\times10^{-16}$ | E001 |
| Bound slack at $p=20$ | median $\approx4\times$, max $5.8\times$, 840 runs | E003 |
| Depth accumulation ratio ($p{=}30/20$) | 1.505 pooled (per-instance $1.54\pm0.26$) vs predicted 1.5 | E011, E033 |
| Small-angle exponents ($\beta$, $\gamma$) | 0.98, 1.92 | E007 |
| Density law | $\lambda/(\beta\gamma^2 m)$ const $\pm40\%$ | E007 |
| $p=30$ small-ramp overlaps, $n=16$--$20$ | 0.71--0.81 (0.71 at $n=20$) | E011 |
| Sublinear-in-$m$ deep leakage | +13--30% vs $m$ +30--58% | E011 |
| Trajectory effective dimension | 2--4 (99% energy) | E009 |
| Exact-compression regret, $p=1$ | 0.028--0.038 ($n=12$--$14$), 0.032--0.055 ($n=16$--$18$) | E002, E010 |
| Regret at $p=3$ ramps | 0.06--0.11 | E002, E010 |
| Regret growth $n=12\to18$, $p=1$ | +0.005 (dense ER) to +0.018 (3-regular), monotone | E002, E010 |
| Leakage-regret family ranking, $p=1$ | $\rho=0.96$ CI $[0.64,0.96]$ ($n{=}12$); $0.86$ $[0.61,0.93]$ ($n{=}14$); pre-committed $\rho\ge0.8$ | E004, E033 |
| Argmax displacement vs regret | pooled $\rho=0.74$ (CI 0.64--0.82) $p=1$; 0.65 (0.54--0.74) $p=3$; robust 0.56--0.74 under all controls | E014, E033 |
| Fidelity confounded | pooled $0.39\to0.07$ within cells ($p1$); pooled $-0.02$ hiding within-cell 0.46 ($p3$); flat-peak $\lvert\rho\rvert\le0.24$ | E014, E033 |
| $n=16$ stability | displacement $\rho$ 0.84/0.76 (normalized conv.; 0.25/0.47 unnormalized); anatomy: 69/70 overprediction, curse 4.2/9.4 | E040 |
| Analytical model off-ER regret, $p=1$ | $\approx$0.01--0.02 (paper quotes 0.005--0.02) | E004 |
| Dense-ER artifact | $\langle C\rangle=93$ on 38 edges, norm $\times7.5$ | E004 |
| Filter backfire | regret $\approx0.01\to0.03$--$0.31$ | E005, E006 |
| Triangle fit paradox | regret $0.054\to0.211$, worse on 136/140 | E012 |
| Gaussian fit inert | MSE $10\times$ better, argmax moves 0/140 | E012 |
| No scalar norm predicts regret | MSE $\rho\approx0.03$ pooled/0.004 demeaned; only displacement $\approx0.7$ | E012, E033 |
| Direction beats size | weighted-norm $\rho=0.74$ vs MSE 0.38; invisible directions free at $\varepsilon=0.5$, aligned harmful at 0.01 | E017 |
| MSE misorders shapes | 105/150 (worse than chance) | E017 |
| Shape saturation in weighted norm | relative error $\ge0.92$; refit doesn't rescue | E017 |
| Contractivity in the wild | exact-$N$ state weight $\le0.984$ on 150/150 | E017 |
| Instrument not recipe | normalizing the analytical model: regret $0.05\to0.16$ | E017 |
| $V_2$ lemmas verified | $10^{-10}$, 280 instances; cubic law $<0.3\%$ | E015 |
| Conditioning correction tracks $\tau^2/m$ | Pearson 0.994 | E015 |
| Quadratic conditioning bound | captures 91--100% (mean 97%) vs 18--64% linear-only; ER limit 0.375 vs 0.352 measured | E035 |
| Weighted-MaxCut extension | identities to $10^{-9}$; capture 92.9--99.6%; cubic law 0.5% ($0.998\pm0.001$); continuous weights collapse classes | E037 |
| Regret certificate | holds 280/280; median tightness 10.5 ($p1$) / 11.3 ($p3$), range 6.4--46 | E036 |
| Signed regret anatomy | $e<0$ at 557/560 evaluations (277/280 pairs both); curse 3.8--4.7$\times$; cancellation 0.6; margin $\approx$40% of remainder | E038 |
| Norm-loss correction (negative) | per-instance $\rho$ $-0.3$ to 0.6; pooled correction better on 0/126, worse on 112 ($p1$) | E039 |
| Depth normalization rule | $p=3$: $0.080\to0.044$ (134/140); $p=1$: 0.031 vs 0.047; $n=16$: $0.090\to0.060$ (32/35) | E036, E040 |
| Sampled-$N$ ($S=10$) match | within 0.011 ($p1$) / 0.016 ($p3$) AR; regret lower by 0.003--0.006, 26/28 cells | E010, E033 |
| Sampled leakage ($S=5$) | median 3.2% relative error | E008 |
| Ceiling validation | $p=1$ gap $\le0.0007$; ramp-grid gap 0.002 (max 0.005); ramp restriction +0.007 (max 0.05) | E034 |
| Transfer, single $n=10$ source | pooled regret 0.007 (median 0.003) vs proxy pooled 0.061 | E018 |
| Transfer, 10-instance mean angles | regret $\approx$0.014 ($p1$) / 0.006 ($p3$); beats proxy 134/140 and 140/140 | E034 |
| Harsh depth ($p=10,20$) | proxy variants 0.07--0.37 (exact 0.07--0.14, sampled 0.11--0.21, analytical 0.18--0.37); transfer $\le0.03$, best-of-3 $\le0.007$ | E019 |
| Compression argmax drift at depth | $\sim$0.1 AR at fidelity 0.7--0.8 | E019, E011 |
| Beyond the wall ($n=22,26$) | transfer wins every cell; analytical improves with $n$ on sparse regular ($0.019\to0.017$) | E019 |
| Analytical P load-bearing | empirical P degrades 3-regular $p=3$: $0.05\to0.27$ | E019 |
| Binned weighted proxy | $\lambda^2\approx\lambda_{\mathrm{struct}}^2+O(1/K)$; $K=64$--$128$: 0.038 vs $\sim$0.032 ($p1$ ER), 0.103 vs 0.105 ($p3$ 3-reg); $O(K^2n)$ | E020 |
| Mild weights | transfer regret up 3--10$\times$ but wins every cell | E021, E022 |
| Scale audit (retraction) | Pareto $G(16,0.5)$ $p1$: proxy 0.012; mean-rescaled transfer 0.090; cost-width-rescaled 0.015 / best-of-3 0.002 | E032 (vs E022--E026) |
| Corner is a grid artifact | continuous refinement residual $<0.006$ from every start; ceiling itself sat 0.020--0.032 low | E029 |
| Guarded surrogate polish | ramp polish +0.05--0.12; full-$2p$ overfits (Pareto $\to0.829$); guard never hurts | E030, E031 |
| Max-3-XOR replication | $p=1$ regret 0.031--0.040; sampled $\ge$ exact 7/8 cells; leakage $0.045\to0.092$ with density; large-angle 0.84--0.95 | E027 |
| Timing wall | exact $N$ 21.8 s at $n=20$; sampled 0.44 s; sweep 0.38 s; GPU ceiling 0.48 s; statevector 0.3 ms | E013 |
| sampled+norm composes (Sec. 6.1) | $p1$ 0.021 vs exact-norm 0.047 (123/7); $p3$ 0.024 vs 0.044 (112/23) | E041 |
| noise is the ingredient (Sec. 6.1) | replicate-averaged loses gain; $S=3$ best (0.012/0.016); converges to exact by $S=300$ | E042 |
| structured noise required at depth (Sec. 6.1) | iid noise wins at $p1$ ($\sigma=0.1$: 0.011) but only hurts at $p3$ (0.036--0.126) | E043 |
| composed recipe replicates at n=16 (Sec. 6.1) | 34/35 both depths; 0.018/0.024 vs exact-norm 0.061/0.060; $S=3$: 0.0067 | E044 |

## If you remember five things

1. **Proxy = QAOA + a blur (projection) after every layer** (Thm 1). Not an
   approximation by definition; an exact compression.
2. **Two errors, independent: compression error (the blur) and model error
   (a wrong $N$).** Most past confusion is mixing them up.
3. **Leakage** $\lambda_\ell$ is the cheap ($O(2^n)$) ruler for compression
   error; it bounds the true-proxy distance (Thm 2), equals a
   *mixer-weighted within-class variance* (Thm 3), and reaches regret
   through the two-point certificate (loose for a measured, interesting
   reason: overprediction + winner's curse + margin).
4. **At small angles, compression error is third order ($\beta\gamma^2$)**
   because the phase separator never leaks and the mixer's first-order leak
   cancels by a MaxCut 2-locality identity. Density sets the constant, and
   the conditioning bound (codegrees, triangles, 4-cycles) pins it to
   within a few percent.
5. **Calibrated honestly, the proxy is not a faster parameter setter.**
   Parameter setting is argmax transfer, which is why MSE fitting
   backfires, why no scalar norm certifies a model, and why plain
   within-family transfer (with a robust cost scale under weights) matched
   or beat every proxy variant at every depth tested. What survives is the
   instrument: sampled $N$ as the best proxy variant at $p\le3$, quantile
   bins for continuous costs, guarded surrogate polishing, the analytical
   argmax beyond simulability, and above all the leakage certificate.
