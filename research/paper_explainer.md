# A plain-language companion to the paper

*(Note: this file is the companion to the **paper draft** specifically. The
living walkthrough of all results, including experiments not yet folded into
the paper, is [explainer.md](explainer.md).)*

*Written for Spencer, 2026-06-13; updated 2026-07-22 after E014/E015/E016 were
folded into the draft and every number was re-audited against the experiment
records. Goal: explain the paper's ideas in order, intuition first and symbols
second, completely enough that you could rewrite the paper from scratch from
this document alone. Nothing here is new science; it is the paper's content,
unpacked. Section 9 gives the paper's skeleton and section 10 the key-numbers
table with experiment provenance.*

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
$Q_\ell(v)$, one complex number per cost $v$, evolved by a recursion whose
coefficients are the **cost-and-distance distribution** $N(v';d,v)$ —
"starting from a typical bitstring of cost $v'$, how many bitstrings of cost
$v$ sit at Hamming distance $d$?"

For the rewrite you need the recursion itself. With mixer matrix elements
$f_d(\beta)=(\cos\beta)^{n-d}(-i\sin\beta)^d$ depending only on Hamming
distance $d$:

$$Q_\ell(v') = \sum_{d,v} f_d(\beta_\ell)\, e^{-i\gamma_\ell v/2}\,
Q_{\ell-1}(v)\, N(v';d,v), \qquad Q_0 \equiv 2^{-n/2},$$

and the empirical $N$ is the class average
$N(v';d,v)=\frac{1}{M_{v'}}\sum_{x\in S_{v'}} n(x;d,v)$. After $p$ layers the
predicted objective is $\langle C\rangle \approx \sum_{v} M_v\,|Q_p(v)|^2\,v$
(one term per attained cost; the paper never needs the $2^n\times2^n$
picture). Each layer costs $O(nm^2)$.

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
Experimentally the bound's slack at $p=20$ is ${\approx}4\times$ in the median
and never above $5.8\times$, and $\sum\lambda_\ell$ is a near-functional
*predictor* of the actual error, not just an upper bound.

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

### 4e. Opening the $V_2$ black box — the two codegree lemmas (E015, now in §4)

The Corollary leaves $V_2$ as an opaque variance. Two exact, machine-verified
lemmas (to $10^{-10}$ on 280 instances) reduce it to pure graph structure:

- **Variance reduction.** Write $\delta_i(y) = c(y\oplus e_i) - c(y)$ for the
  single-flip cost changes and $T(y) = \sum_i \delta_i(y)^2$. Then
  $s_2 = (n-8)c^2 + 4mc + T$, and the first two terms are functions of the
  cost alone, so they vanish inside a within-class variance:
  **$V_2$ is just the within-class variance of $T$** — only *squared
  single-flip cost changes* matter. (Proof: expand $(c+\delta_i)^2$, apply the
  neighbor-sum Lemma to $\sum_i\delta_i = 2m-4c$.)
- **Codegree form.** In spin variables $\sigma_j=(-1)^{y_j}$, one line of
  algebra gives $\delta_i = \sigma_i\sum_{j\sim i}\sigma_j$, hence
  $T = \sum_{j,k} A_{jk}\sigma_j\sigma_k$ where $A_{jk}$ counts **common
  neighbors** of vertices $j,k$ (the codegree). Over uniform bitstrings this
  gives *exactly* $\mathbb{E}[T]=2m$ and
  $\mathrm{Var}(T) = 2\sum_{j\ne k}A_{jk}^2$ — a closed form in $O(n^2+m)$,
  no $2^n$ enumeration. For $G(n,p)$ it evaluates to $\to 8p^2m^2$.

**The conditioning story (the mechanism behind the density law).** $V_2$ is
the *cost-conditioned* version of $\mathrm{Var}(T)$. Unconditionally,
$\sqrt{\mathrm{Var}(T)}/m$ grows with density $p$. But conditioning on the
cost removes a fraction that is **triangle-driven** (only codegree pairs that
are also edges, i.e. triangles, feel the conditioning) and that fraction also
grows with density. The two nearly cancel, leaving the measured flat law
$\sqrt{V_2}\propto m$. Empirically the removed piece
$\mathrm{Var}(T)-V_2$ tracks $\tau^2/m$ ($\tau$ = triangle count) at Pearson
$0.994$. That correction is no longer an open prefactor: the quadratic
projection bound of §4f pins it rigorously (it captures 90–100%, mean 97%,
of the exact correction, gives an explicit $G(n,p)$ limit formula, and
leaves at most a 10% residual). The operative identity is
$\mathrm{Cov}(T,S^2)=\mathrm{Var}(T)$, also recorded in §7 item 2; the
triangle channel alone is only part of the story, which is exactly why the
crude linear estimate undershot. Slogan for the paper: *"density drives
compression error" sharpens to "squared codegrees, conditioned on cost, drive
it."*

### 4f. The quadratic conditioning bound (Proposition, E019)

This is the piece that turns the density law from an observation into a
formula. Everything runs through the **edge sum** $S = m - 2c$, a linear
function of the cost, so conditioning on $c$ is the same as conditioning on
$S$. Four exact identities over uniform bitstrings (with $\tau$ the triangle
count and $c_4$ the 4-cycle count):

$$\mathbb{E}[S^3]=6\tau,\quad \mathrm{Var}(S^2)=2m^2-2m+24c_4,\quad
\mathrm{Cov}(T,S)=6\tau,\quad \mathrm{Cov}(T,S^2)=\mathrm{Var}(T).$$

The last one is the operative one: it is a **quadratic channel present even
in triangle-free graphs**, which is why the crude linear estimate
($36\tau^2/m$) undershot the correction so badly (it captures only 18–64%).

**The Gram bound.** Project $T$ onto $\mathrm{span}\{S,\,S^2-m\}$. The
variance of that $L^2$ projection is $v^\top G^{-1} v$ with
$v = (6\tau,\ \mathrm{Var}(T))$ and $G$ the Gram matrix built from
$\mathrm{Var}(S)=m$, $\mathbb{E}[S^3]=6\tau$, and $\mathrm{Var}(S^2)$. Since
any projection onto a subspace is dominated by the projection onto *all*
functions of $S$, which is exactly $\mathrm{Var}(\mathbb{E}[T\mid S]) =
\mathrm{Var}(\mathbb{E}[T\mid c])$, we get

$$\mathrm{Var}\bigl(\mathbb{E}[T\mid c]\bigr) \ge v^\top G^{-1} v,
\qquad\text{equivalently}\qquad
V_2 \le \mathrm{Var}(T) - v^\top G^{-1} v.$$

**Proof idea.** Expand $T - 2m = 2\sum_{j<k}A_{jk}\sigma_j\sigma_k$ (the
codegree form) and $S = \sum_{(a,b)\in E}\sigma_a\sigma_b$ in spin
variables; over uniform spins **only even monomials survive**, so every
moment is a graph count. $\mathbb{E}[\sigma_j\sigma_k S]$ is an edge
indicator, and summing $A_{jk}$ over edges counts each **triangle at each
of its three edges**, twice per orientation, giving $6\tau$.
$\mathbb{E}[\sigma_j\sigma_k S^2] = 2A_{jk}$ (**ordered two-edge paths**
from $j$ to $k$), and the $A_{jk}$-weighted sum of that is exactly
$4\sum_{j<k}A_{jk}^2 = \mathrm{Var}(T)$. $\mathbb{E}[S^3]$ counts ordered
triangles; $\mathbb{E}[S^4]$ counts edge pairings plus oriented 4-cycles.

**Why it matters practically:** every ingredient ($m$, $\tau$, $c_4$,
codegrees) is computable from the graph in **polynomial time**; no $2^n$
enumeration anywhere. Across 140 instances of the seven families the bound
captures 90–100% (mean 97%) of the exact conditioning correction. For
$G(n,p)$ it evaluates asymptotically to

$$V_2/\mathrm{Var}(T) \;\to\; 1 - \frac{p(1+4p-2p^2)}{1+6p^2-4p^3},$$

giving $0.375$ at $p=\tfrac12$ against $0.352$ measured at $n=12$–$14$. The
empirical density law becomes an explicit formula, up to a residual of at
most 10%.

### 4g. Weighted MaxCut: the theory travels, the object doesn't

All of §4c–4f extends verbatim to integer-weighted MaxCut under mechanical
substitutions (Remark "Weighted MaxCut" in the paper, E021):

- **Neighbor-sum Lemma:** $m \to W$ (total weight).
- **Second moments:** replace $m$ by $\sum_e w_e^2$ and pick up a
  $\sum_e w_e^4$ term: $\mathbb{E}[T] = 2\sum_e w_e^2$,
  $\mathrm{Var}(S) = \sum_e w_e^2$, and
  $\mathrm{Var}(S^2) = 2\bigl(\sum_e w_e^2\bigr)^2 - 2\sum_e w_e^4 + 24c_4^w$.
- **Codegrees:** $A_{jk} \to \sum_i w_{ij}w_{ik}$.
- **Triangle and 4-cycle counts** $\to$ the corresponding weight products.

All identities were machine re-verified to $10^{-9}$ on 140 weighted
instances (bound capture again 93–99.6%, cubic law to 0.5%). The caveat that
delimits the whole enterprise: generic **continuous** weights collapse every
cost class to a bitstring-complement pair, so the compression itself becomes
vacuous. Cost degeneracy is a requirement of the *object*, not of the
theory.

### 4h. From state error to parameter error: the two-point regret certificate (E020)

Everything so far bounds the **state** error, but parameter setting consumes
an **argmax**. The paper bridges the two with a certificate evaluated at
just two points (the subsection "From state error to parameter error",
Proposition "Two-point regret certificate").

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
the two relevant points, and both sit in the low-leakage corner in practice.

**Proof idea (three moves).** (1) Decompose the regret into three terms,
$[F-\hat F](\theta^*) + [\hat F(\theta^*)-\hat F(\hat\theta)] +
[\hat F-F](\hat\theta)$; the **middle term is nonpositive** because
$\hat\theta$ maximizes $\hat F$. (2) For the two survivors, use the identity
(valid for normalized states and any scalar $a$)
$\langle\psi|C|\psi\rangle - \langle\chi|C|\chi\rangle =
\langle\psi-\chi|(C-a)\psi\rangle + \langle\chi|(C-a)(\psi-\chi)\rangle$;
**center at $a=\langle C\rangle_\psi$** and apply Cauchy–Schwarz, using
$\lVert(C-a)\chi\rVert^2 = \sigma_\chi^2+\delta^2$. (3) Bound the distance
via **Theorem 2**: $\lVert\psi-\phi\rVert \le \sum_\ell\lambda_\ell$, and
normalizing $\phi$ costs at most another
$1-\lVert\phi\rVert \le \lVert\psi-\phi\rVert$, hence the factor 2.

**The unnormalized variant.** The certificate as stated covers the
normalized proxy objective, which is also the convention recommended at
depth (§5). For the unnormalized objective preferred at $p=1$ it acquires
one extra term, $|a|\,(1-\lVert\phi\rVert^2) = |a|\sum_\ell\lambda_\ell^2$
with $a = \langle C\rangle_\psi$.

**How to read it.** The certificate holds on all 280 instance-depth pairs
tested (machine-checked) but is honestly loose: median tightness is
$10$–$11\times$ over the actual regret, informative only where accumulated
leakage is a few percent. And that looseness *is* the argmax-transfer story
of §5, quantified: regret is small not because the proxy landscape is
pointwise accurate, but because its errors at $\theta^*$ and $\hat\theta$
nearly cancel, a cancellation a triangle inequality cannot see. Sharpening
it (a correlated-error bound) is the natural next theory question.

### 4d. The "why homogeneity happens" ladder

The paper's closing §4 subsection assembles the above into a narrative
answer to "why do equal-cost bitstrings end up with equal amplitudes?":

1. The initial state $\ket{+}^{\otimes n}$ is **exactly** homogeneous.
2. The phase separator preserves homogeneity **exactly** (Step 1).
3. The mixer is the only source of leakage, and Theorem 3 says only
   *within-class variance* of the weighted neighborhood enters.
4. At small $\beta$ that variance is dominated by single flips, where the Lemma
   makes the first moment an exact function of cost — **leading-order protection
   on any graph** (it's 2-locality, not randomness).
5. Randomness buys only the *next* rung: $V_2$ self-averages over the
   exponentially large cost classes for random-like instances. **This is the
   one rung that is measured, not proved** — though §4e–4f now pin the
   density law to within a few percent, leaving only the small residual and
   a rigorous self-averaging statement.
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
→ **At $p=1$, yes, and the ranking holds:** regret $0.028$–$0.055$ at $p=1$
(roughly doubling at $p=3$ ramps to $0.06$–$0.11$), positive value-added on
every family, mild growth with $n$ ($+0.003$–$0.017$ from $n=14$ to $18$).
Ranking families by leakage *at the proxy's own chosen angles* reproduces
their ranking by regret (Spearman $\rho=0.96$ at $n{=}12$, $0.86$ at
$n{=}14$). **Honest limit:** the ranking is a $p=1$ statement over seven rank
points; along proxy-chosen $p=3$ ramps it collapses to $\rho\approx0$
because those schedules equalize leakage across families.

**Q: At depth, is regret a fidelity problem at all? (E014 — new)**
→ **No — it is an argmax-transfer problem, and fidelity decouples.**
Recomputing true-landscape geometry over 140 instances at both depths:
regret correlates with the proxy's **argmax displacement** from the true
optimum at both depths (pooled Spearman $\rho=0.74$, 95% CI $0.64$–$0.82$,
at $p=1$; $0.65$, CI $0.54$–$0.74$, at $p=3$; and robustly
$\rho=0.56$–$0.74$ under every control: family-and-size demeaning,
within-cell averaging, cluster bootstrap). The state-fidelity deficit is
**family-confounded both ways**: its pooled $\rho=0.39$ at $p=1$ drops to
$0.07$ within family-size cells, while the $p=3$ pooled value of $-0.02$
hides a moderate within-cell $0.46$ (proxy-chosen schedules compress
fidelity into a narrow band across families). Landscape flat-peak
robustness predicts nothing ($|\rho|\le0.24$ under all controls). So the
leakage calculus bounds the *state* error while parameter quality lives in
*parameter space* — which is exactly why a lower-fidelity model with a
better-placed argmax (the analytical $N$ off-ER) beats the exact
compression. This is the result that lets the paper narrate depth honestly.

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

**Q: how does the proxy compare to plain parameter transfer? (E018 — new)**
→ **Transfer wins on these ensembles.** The grid ceilings themselves are
validated (continuous refinement gains at most 0.0007 at $p=1$ and ~0.002 on
the $p=3$ ramp grid; dropping the linear-ramp restriction gains another
~0.007). But taking the mean true-optimal angles of ten brute-forced
ER(0.5) $n=12$ instances and applying them verbatim to every instance,
family, and size beats the exact-compression proxy on 134/140 instances at
$p=1$ and 140/140 at $p=3$. Parameter concentration is that strong on random
ensembles. The honest conclusion (now in the paper): instance-specific
parameter setting pays only where concentration fails (structured,
heterogeneous instances) or where nothing can be simulated (then only the
analytical $N$ exists anyway).

**Q: why normalize the proxy objective at depth? (E020 — new)**
→ **Because the unnormalized objective carries a leakage bias.** Along the
depth grid the compressed norm varies substantially from schedule to
schedule, so the unnormalized objective **conflates "high predicted value"
with "low leakage"** and systematically favors small-angle schedules.
Dividing by the tracked compressed norm (free, from Theorem 2's norm
bookkeeping) removes the bias exactly where leakage accumulates. Effect:
pooled $p=3$ regret drops $0.080\to0.044$ (better on 134/140 instances); at
$p=1$ the unnormalized convention keeps a mild edge ($0.031$ vs $0.047$),
so the rule is: **normalize at depth, keep unnormalized at $p=1$.**

**Q: so what should a practitioner actually do?**
→ **Follow the three-branch decision rule** (stated verbatim in the paper's
"How to set parameters" paragraph at the end of §5.4):

1. **If the ensemble concentrates and a solved source instance exists** (the
   random families tested here, at simulable sizes): **transfer** its
   angles. Transfer beat every per-instance method tested (E018).
2. **Otherwise, where statevector passes are affordable:** the recipe is
   **sampled $N$ ($S\approx10$ bitstrings per cost class) + the empirical
   cost distribution, with the normalized objective at depth.** It matches
   exact parameter setting within $0.02$ AR, regret no worse,
   family-agnostic, no dense-graph artifact. **Honest scope:** its cost is
   $O(S\,m\,2^n)$, still exponential, useful only where statevector passes
   are affordable but the $O(4^n)$ exact distribution is not (roughly
   $n\lesssim30$). Whether $N$ has a polynomial Monte-Carlo estimator is
   open.
3. **Beyond classical reach, the analytical $N$ is the only option:** trust
   its argmax, never its values, and distrust it on dense Erdős–Rényi
   graphs.

**Q: stepping back, when does the proxy work at all?**
→ This is the three-part answer that now *opens* the paper's Discussion
("So when does it work?"):

1. The **compression** is faithful whenever density-controlled leakage is
   small at the schedules worth running, which the argmax's self-selection
   makes generic (the §4 ladder).
2. Its **parameter setting** works whenever the argmax transfers, which held
   on every tested family and failed only for the analytical $N$ on dense
   Erdős–Rényi.
3. It **adds value over the cheapest alternative** (transferring angles from
   one solved instance) only where concentration fails or no solved source
   exists: structured, heterogeneous, or beyond-classical instances.

---

## 6. What is proved vs. measured (so you can defend it)

- **Proved (exact, machine-verified):** Theorem 1 (proxy = compression),
  Theorem 2 (telescoping bound + norm identity), Theorem 3 (variance identity),
  the Lemma (neighbor-cost sum), the Corollary's cancellation (first-order leak
  is identically zero), and the two $V_2$ lemmas (variance reduction; codegree
  form, including the exact $G(n,p)$ expectation). These hold for *any* graph;
  no randomness needed.
- **Measured, not proved:** that $V_2$ self-averages so that leakage tracks
  density across random-graph families; the family rankings; all regret numbers;
  the sampled-$N$ recipe's quality. These are empirical over the tested range
  ($n\le20$, $p\le30$, 7 families, 20–30 instances/cell, fixed seeds).
- **The honest seams** (the paper states all of these):
  - The leakage→regret ranking is a **7-point Spearman** that is only resolvable
    at $p=1$/small $n$; it fades at scale. Don't oversell it.
  - Experiments are **unweighted MaxCut**, $n\le20$; the theory itself is
    weighted-ready (E021), but continuous weights dissolve the compression
    (cost classes collapse to complement pairs).
  - Trajectory PCA is a *diagnostic*, not a method (it needs the states the proxy
    exists to avoid computing).
  - The recipe's cost is still exponential.

---

## 7. The open doors — one now measured shut, one narrowed

1. **Instance-adapted low-rank frames — the cheap version fails (E016, now in
   the Discussion).** The trajectory really is ~4-dimensional (an oracle
   rank-4 PCA frame captures ~0.99 of a depth-20 trajectory), so a good
   instance-adapted frame would easily beat the $(m{+}1)$-dim cost-class
   frame. But a rank-4 frame built from the first five layers is nearly
   **orthogonal** to the trajectory's true late-time subspace (largest
   principal angle $71°$–$88°$, growing with ramp size) and captures *less*
   than the zero-cost cost-class frame. The low-rank subspace **rotates**
   across depth, so it cannot be discovered from a cheap prefix. Any frame
   that beats the compression must model that rotation explicitly. This is a
   measured obstacle, reported as future-work guidance, not a method.
2. **The $V_2$ density law — now pinned (E019, new Proposition in §4).**
   The missing piece was found: projecting $T$ onto $\{S, S^2\}$ of the edge
   sum $S=m-2c$ gives a rigorous, polynomial-time bound
   $V_2 \le \mathrm{Var}(T) - v^\top G^{-1} v$ built from four exact
   identities; the crucial one is $\mathrm{Cov}(T, S^2)=\mathrm{Var}(T)$, a
   quadratic channel present even in triangle-free graphs, which is exactly
   why the crude linear estimate ($36\tau^2/m$) undershot. The bound
   captures 90-100% (mean 97%) of the exact correction on 140 instances,
   and for $G(n,p)$ evaluates to
   $V_2/\mathrm{Var}(T)\to 1-p(1+4p-2p^2)/(1+6p^2-4p^3)$ (0.375 at $p=1/2$
   vs. 0.352 measured). Remaining open: only the few-percent residual and a
   rigorous self-averaging statement.
3. **A polynomial-time estimator of $N$** (or of leakage) would extend the
   recipe beyond the classically simulable regime; nothing rules it out.
4. **The theory-to-regret bridge exists but is loose (E020, new §4
   Proposition).** A two-point certificate
   (regret $\le \varepsilon(\theta^*)+\varepsilon(\hat\theta)$, each
   $\varepsilon$ controlled by leakage and cost variances) is proved and
   machine-verified; median tightness is only 10-11x, and that looseness IS
   the argmax-transfer story quantified: the proxy's landscape errors at the
   two relevant points nearly cancel, which a triangle inequality cannot
   see. Sharpening it (a correlated-error bound) is the natural next theory
   question. Bonus practical rule discovered on the way: **normalize the
   proxy objective at depth** (divide by the tracked compressed norm; halves
   p=3 regret at zero cost; keep unnormalized at p=1).

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
| $\delta_i(y)$ | "single-flip change" | $c(y\oplus e_i)-c(y)$; spin form $\sigma_i\sum_{j\sim i}\sigma_j$ |
| $T(y)$ | "flip-energy sum" | $\sum_i \delta_i^2$; $V_2$ = its within-class variance |
| $A_{jk}$ | "codegree" | # common neighbors of $j,k$; $\mathrm{Var}(T)=2\sum_{j\ne k}A_{jk}^2$ |
| $\tau$ | "triangles" | drives the conditioning correction $\mathrm{Var}(T)-V_2 \propto \tau^2/m$ |

---

## 9. The paper's skeleton (for rewriting it from scratch)

The draft (`qce2027_paper.tex`, ~12 pages single-column) is organized so that
each section owns one move. If you rewrite, this is the load-bearing order:

**The abstract** (~200 words) is shaped thesis-first, no suspense: (i) the
thesis in the first two sentences (the proxy is an *exact orthogonal
compression*, not a substitution heuristic); (ii) the error calculus it
yields (leakage, the telescoping bound with its median $4\times$ slack, the
variance identity, and the two-point certificate transferring state error to
regret); (iii) the density law, derived within ten percent from exact
codegree and cycle identities; (iv) four experimental headlines (regret
tracks argmax displacement and no distribution- or amplitude-space norm;
sampled $N$ matches exact within $0.02$ AR; parameter transfer beats every
per-instance method on concentrated ensembles; normalizing the objective
nearly halves depth regret at no cost).

1. **Introduction.** The hook is the two open questions (what is the proxy /
   when does it work), quoted against Sud et al.'s own "no longer unitary,
   analogues of amplitudes" language. Contributions list: exactness, error
   calculus, small-angle structure (with the $V_2$ lemmas), experimental
   anatomy, practical findings. A glossary table carries the paper-specific
   terminology.
2. **Background and related work.** One short QAOA subsection; one proxy
   subsection (the recursion, coefficients $N$); related work does the
   novelty positioning: exact symmetry reductions (Shaydulin, Tsvelikhovskiy)
   as the exact special case, lumpability/aggregated Markov chains and model
   order reduction (Buchholz, Antoulas) as the classical mathematics never
   before connected to the proxy, pseudo-Boltzmann $p=1$ states (Diez-Valle)
   and mean-field AOA (Misra-Spieldenner) as the other surrogate families,
   parameter transfer (Brandao, Galda, Sureshbabu, utility-scale 2026) as the
   competing strategy, and Kruger-Mauerer as the closest landscape-surrogate
   neighbor (they approximate the landscape; we ask which error norms certify
   argmax transfer). Three references deliberately live *not* here but in
   §5.4 (khairy2020, shaffer2023, grips2024): the fitted-shape paradox
   re-tests the collaboration's documented G-RIPS 2024 proposal, framed as
   the surrogate-fit pattern of khairy2020/shaffer2023 transplanted into
   distribution space.
3. **Section 3: exactness.** Setup ($S_v$, $M_v$, class states, $P$), the
   empirical $N$, Theorem 1 (one proxy step = layer + projection; proof is a
   two-line induction: group the mixer sum by distance and cost, class-average
   turns $n$ into $N$), the transfer-matrix form
   $T = D^{-1/2}(PUP)|_{\mathcal{H}_{hom}}D^{1/2}$, three delimiting remarks
   (same-instance $N$ only; unattained costs inert; do not renormalize), and
   the norm certificate (Prop. 1: norm inflation proves model error).
4. **Section 4: the error calculus.** Theorem 2 (telescoping + Pythagoras),
   the "leakage is an instrument, $O(2^n)$" paragraph, Theorem 3 (variance
   identity + the wrong-norm warning it implies), then MaxCut small-angle
   structure: neighbor-sum Lemma, cubic Corollary, the two $V_2$ lemmas with
   proofs, the conditioning remark, the Proposition "Quadratic conditioning
   bound" (`prop:quadcond`; §4f here), the Remark "Weighted MaxCut"
   (`rem:weighted`; §4g here), the subsection "From state error to parameter
   error" with the Proposition "Two-point regret certificate"
   (`prop:certreg`; §4h here), and the "why homogeneity happens" ladder.
5. **Section 5: experimental anatomy.** Setup (7 families, $n=12$–$20$,
   grids, regret/value-added/ceiling definitions, SE discipline; now also
   the ceiling validation (Nelder–Mead refinement gains $\le0.0007$ AR at
   $p=1$ and $\approx0.002$ on the $p=3$ ramp grid, releasing the ramp
   restriction a further $\approx0.007$) and the normalization convention
   (argmaxes use the unnormalized objective unless stated; normalized drops
   pooled $p=3$ regret $0.080\to0.044$)). 5.2:
   leakage maps, density law, bound tightness, depth scaling, trajectory PCA.
   5.3: regret table, mild-$n$ growth, the $p=1$ leakage-regret ranking
   (honestly scoped), ending with the transfer-baseline paragraph (E018).
   5.4: model error (analytical robustness + dense-ER
   artifact, filter negatives, fitted-shape paradox table, E014
   argmax-vs-fidelity, the normalization-rule paragraph, sampled-$N$ recipe,
   the "How to set parameters" decision rule, timing table).
6. **Discussion.** Opens with "So when does it work?" (the three-part answer,
   see the last Q of §5); then what the compression view buys; limits; E016
   (the rotating subspace kills cheap instance-adapted frames); the $V_2$
   open piece is now only the small residual plus a rigorous self-averaging
   statement.
7. **Code and data availability.** Everything traces to
   `research/experiments/E001–E021`; claim-to-experiment mapping is in LaTeX
   comments (`% E00x` next to each claim).

Style rules the draft follows: every quantitative sentence carries a `% E00x`
comment; negative results are reported as findings, not buried; each
theorem's scope (any graph vs. random-like) is stated where it is used.

## 10. Key numbers, with provenance (post-audit, 2026-07-22)

| Claim | Number | Experiment |
|---|---|---|
| Thm 1 exactness verified | $2.5\times10^{-16}$ | E001 |
| Bound slack at $p=20$ | median $\approx4\times$, max $5.8\times$, 840 runs | E003 |
| Depth accumulation ratio ($p{=}30/20$) | $1.505$ vs. predicted $1.5$ | E011 |
| Small-angle exponents ($\beta$, $\gamma$) | $0.98$, $1.92$ | E007 |
| Density law | $\lambda/(\beta\gamma^2 m)$ const $\pm40\%$ | E007 |
| $p=30$ small-ramp overlaps, $n=16$–$20$ | $0.81/0.76/0.71$ | E011 |
| Sublinear-in-$m$ deep leakage | $+11$–$33\%$ vs. $m$ $+25$–$57\%$ | E011 |
| Trajectory effective dimension | $2$–$4$ (99% energy; ~13 at large ramps) | E009 |
| Exact-compression regret, $p=1$ | $0.028$–$0.038$ ($n\le14$), $0.032$–$0.055$ ($n=16$–$18$) | E002, E010 |
| Regret at $p=3$ ramps | $0.063$–$0.107$ | E002, E010 |
| Regret drift $n=14\to18$ | $+0.003$–$0.017$ | E010 |
| Leakage-regret family ranking, $p=1$ | Spearman $0.96$ ($n{=}12$), $0.86$ ($n{=}14$); $\approx0$ at $p=3$ | E004 |
| Analytical model off-ER regret | $0.01$–$0.02$ | E004 |
| Dense-ER artifact | $\langle C\rangle=93$ on 38 edges, norm $\times7.5$ | E004 |
| Filter backfire | regret $0.01\to0.07$–$0.31$ | E005, E006 |
| Triangle fit paradox | regret $0.054\to0.211$, worse on 136/140 | E012 |
| Gaussian fit inert | median MSE $10\times$ better, argmax moves 0/140 | E012 |
| Norms vs. regret | MSE $\rho\approx-0.1$; amplitude $-0.2$; landscape $+0.1$; argmax displacement $\approx0.7$ | E012 |
| Raw analytical slice sums | mean $\sim10^9$, up to $2\times10^{10}$ (vs. $2^n$) | E012 |
| Argmax vs. fidelity at depth | displacement pooled $\rho=0.74$ (CI $0.64$–$0.82$) at $p=1$, $0.65$ (CI $0.54$–$0.74$) at $p=3$, robust $0.56$–$0.74$ under all controls; fidelity pooled $0.39\to0.07$ within cells ($p{=}1$), $-0.02$ pooled hiding within-cell $0.46$ ($p{=}3$); robustness $\lvert\rho\rvert\le0.24$ | E014, E017 |
| $V_2$ lemmas verified | $10^{-10}$, 280 instances; cubic law $<0.3\%$ | E015 |
| Conditioning correction | $\propto\tau^2/m$, Pearson $0.994$ | E015 |
| Prefix-frame failure | principal angle $71°$–$88°$; oracle captures $0.99$ | E016 |
| Quadratic conditioning bound | captures 90–100% (mean 97%) of $\mathrm{Var}(\mathbb{E}[T\mid c])$; ER limit $0.375$ at $p=1/2$ | E019 |
| Ceiling validation | $p=1$ gap $\le0.0007$; ramp gap $\approx0.002$; ramp restriction $\approx0.007$ | E018 |
| Transfer baseline | beats proxy 134/140 ($p1$), 140/140 ($p3$); regret $\approx0.014/0.008$ | E018 |
| Statistics hardening | argmax $\rho$ 0.56–0.74 under all controls; fidelity family-confounded both ways | E017 |
| Regret certificate | holds on 280/280; median tightness 10–11x | E020 |
| Depth normalization rule | $p=3$ regret 0.080→0.044 (134/140 better); $p=1$: 0.031 vs 0.047 | E020 |
| Weighted MaxCut extension | all identities hold; capture 93–100%; cubic law to 0.4%; continuous weights collapse classes to complement pairs | E021 |
| Sampled-$N$ ($S=10$) match | within $0.011$ ($p1$) / $0.016$ ($p3$) AR; regret lower in 26/28 cells | E010 |
| Sampled leakage ($S=5$) | median $3.2\%$ relative error | E008 |
| Timing wall | exact $N$: $21.8$ s at $n=20$; sampled: $0.44$ s; brute-force ceiling: $0.48$ s | E013 |
| Pipeline crossover | $n\approx22$–$24$ (extrapolated) | E013 |

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
5. **Parameter setting is argmax-transfer, not state-approximation, and the
   decision rule has three branches.** Transfer angles from one solved
   instance wherever the ensemble concentrates and such a source exists
   (transfer won on every tested random family); otherwise use sampled $N$ +
   the empirical cost distribution, with the normalized objective at depth;
   beyond classical reach the analytical $N$ is the only option (trust its
   argmax, never its values, avoid dense ER). Argmax-transfer is also why
   fitting $N$ by MSE backfires and why the analytical model works far
   outside its derivation domain.
