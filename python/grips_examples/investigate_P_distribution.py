"""
Investigate the effect of using real P(c') vs. binomial P(c') on QAOA proxy
parameter setting.

The QAOA proxy expectation is: <C> = 2^n * sum_c' P(c') * |Q_p(c')|^2 * c'

PaperProxy uses P(c') = Binomial(m, 0.5). This experiment tests whether using
the *real* P(c') (empirical histogram from a specific graph instance) finds
better QAOA parameters.

Key insight: Julia's QAOA_proxy_multi computes Q states from the homodist N,
and expectation() takes P as a separate vector. So we compute Q once, then
evaluate expectations with two different P distributions.

Proxy configurations tested:
  1. Paper N + Paper P  — fully class-level (self-consistent analytical proxy)
  2. Paper N + Real P   — inconsistent hybrid (analytical N, instance P)
  3. Real N  + Real P   — fully instance-specific (requires compute_real_homodist)
  4. Avg N  + Avg P     — batch-averaged instance distributions (requires compute_averaged_distributions)
"""

import numpy as np
import networkx as nx
import qokit.maxcut as mc
from grips.QAOA_proxy_interface import jl
from grips.real_distribution import get_costs
from grips.QAOA_simulator import get_expectation
from grips.solve_maxcut_exact import maxcut

# Bring unexported Julia function into scope
jl.seval("using JuliaQAOA: cpu_compute_homodist")

# ===== Configuration =====
num_nodes = 10
edge_probability = 0.75
seeds = range(5)
num_gamma = 100
num_beta = 50
p = 1  # QAOA depth
compute_real_qaoa = False  # Expensive: full real QAOA grid sweep
compute_real_homodist = False  # Expensive: instance-specific N(c';d,c) + Real P
compute_averaged_distributions = False  # Average N and P across all seeds; requires compute_real_homodist
generate_plots = True  # Set False for text-only output (no matplotlib)

# Build gamma/beta grid once (shared across all seeds)
gammas_1d = np.linspace(0, 2, num_gamma)  # in pi-units (QOKit convention: phase gate is exp(-iγC/2))
betas_1d = np.linspace(0, 0.5, num_beta)  # in pi-units
gamma_mesh, beta_mesh = np.meshgrid(gammas_1d, betas_1d)
gammas_flat = gamma_mesh.flatten().reshape(-1, 1)
betas_flat = beta_mesh.flatten().reshape(-1, 1)
num_pairs = gammas_flat.shape[0]


def run_single_seed(seed):
    """Run the full experiment for one graph seed. Returns a dict of results."""
    print(f"\n{'='*65}")
    print(f"  Seed {seed}")
    print(f"{'='*65}")

    # Step 1: Generate graph and compute real P(c')
    graph = nx.erdos_renyi_graph(num_nodes, edge_probability, seed=seed)
    num_edges = graph.number_of_edges()
    print(f"  n={num_nodes}, m={num_edges}")

    costs = get_costs(graph)
    real_P = np.bincount(np.round(costs).astype(int), minlength=num_edges + 1) / len(costs)
    real_P = real_P[: num_edges + 1]

    # Step 2: Create Julia PaperProxy and compute analytical homodist
    proxy = jl.PaperProxy(num_edges, num_nodes, edge_probability)
    N_paper = jl.cpu_compute_homodist(proxy)

    paper_P = np.array([float(jl.P_cost_distribution(proxy, c)) for c in range(num_edges + 1)])

    # Step 3: Run QAOA_proxy_multi with Paper homodist
    print(f"  Running proxy with Paper N ({num_pairs} grid points)...")
    Qs_paper = jl.QAOA_proxy_multi(N_paper, gammas_flat, betas_flat)
    Q_paper_final = Qs_paper[p]

    # Step 4: Compute proxy expectations (Paper N + Paper P, Paper N + Real P)
    paper_expectations = jl.expectation(Q_paper_final, paper_P, num_nodes).to_numpy().flatten()
    realP_expectations = jl.expectation(Q_paper_final, real_P, num_nodes).to_numpy().flatten()

    paper_landscape = paper_expectations.reshape(num_beta, num_gamma)
    realP_landscape = realP_expectations.reshape(num_beta, num_gamma)

    # Step 5: Instance-specific homodist (optional, expensive: O(2^2n))
    if compute_real_homodist:
        print(f"  Computing instance-specific N(c';d,c) via Julia...")
        N_real = jl.get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_nodes)

        print(f"  Running proxy with Real N ({num_pairs} grid points)...")
        Qs_real = jl.QAOA_proxy_multi(N_real, gammas_flat, betas_flat)
        Q_real_final = Qs_real[p]

        realNP_expectations = jl.expectation(Q_real_final, real_P, num_nodes).to_numpy().flatten()
        realNP_landscape = realNP_expectations.reshape(num_beta, num_gamma)
    else:
        realNP_expectations = None
        realNP_landscape = None

    # Step 6: Real QAOA landscape (optional, expensive)
    ising_model = mc.get_maxcut_terms(graph)
    if compute_real_qaoa:
        print(f"  Computing real QAOA landscape...")
        real_qaoa_expectations = np.zeros(num_pairs)
        for i in range(num_pairs):
            gamma_rad = np.array([gammas_flat[i, 0] * np.pi])
            beta_rad = np.array([betas_flat[i, 0] * np.pi])
            real_qaoa_expectations[i] = get_expectation(num_nodes, ising_model, gamma_rad, beta_rad)
        real_qaoa_landscape = real_qaoa_expectations.reshape(num_beta, num_gamma)
    else:
        real_qaoa_expectations = None
        real_qaoa_landscape = None

    # Step 7: Find optima and evaluate real QAOA at each proxy's best params
    paper_best_idx = np.argmax(paper_expectations)
    realP_best_idx = np.argmax(realP_expectations)

    c_opt, _ = maxcut(graph)

    def real_qaoa_at_idx(idx):
        """Get real QAOA expectation at a grid index (from grid or point eval)."""
        if real_qaoa_expectations is not None:
            return real_qaoa_expectations[idx]
        return get_expectation(
            num_nodes,
            ising_model,
            np.array([gammas_flat[idx, 0] * np.pi]),
            np.array([betas_flat[idx, 0] * np.pi]),
        )

    real_at_paper = real_qaoa_at_idx(paper_best_idx)
    real_at_realP = real_qaoa_at_idx(realP_best_idx)
    apx_paper = real_at_paper / c_opt
    apx_realP = real_at_realP / c_opt

    if compute_real_homodist:
        realNP_best_idx = np.argmax(realNP_expectations)
        real_at_realNP = real_qaoa_at_idx(realNP_best_idx)
        apx_realNP = real_at_realNP / c_opt
    else:
        realNP_best_idx = None
        real_at_realNP = None
        apx_realNP = None

    if compute_real_qaoa:
        real_qaoa_best_idx = np.argmax(real_qaoa_expectations)
        real_at_qaoa = real_qaoa_expectations[real_qaoa_best_idx]
        apx_qaoa = real_at_qaoa / c_opt
    else:
        real_qaoa_best_idx = None
        real_at_qaoa = None
        apx_qaoa = None

    # Print per-seed result
    parts = [f"  c_opt={c_opt}, Paper N+P={apx_paper:.4f}, Paper N+Real P={apx_realP:.4f}"]
    if apx_realNP is not None:
        parts.append(f"Real N+P={apx_realNP:.4f}")
    if apx_qaoa is not None:
        parts.append(f"QAOA opt={apx_qaoa:.4f}")
    print(", ".join(parts))

    # Winner determination
    candidates = {"Paper N+P": apx_paper, "Paper N+Real P": apx_realP}
    if apx_realNP is not None:
        candidates["Real N+P"] = apx_realNP
    best_name = max(candidates, key=candidates.get)
    best_val = candidates[best_name]
    tied = [k for k, v in candidates.items() if v == best_val]
    winner = " / ".join(tied) if len(tied) > 1 else best_name
    print(f"  >> Winner: {winner}")

    return {
        "seed": seed,
        "num_edges": num_edges,
        "c_opt": c_opt,
        "apx_paper": apx_paper,
        "apx_realP": apx_realP,
        "apx_realNP": apx_realNP,
        "apx_qaoa": apx_qaoa,
        "best_gamma_paper": gammas_flat[paper_best_idx, 0],
        "best_beta_paper": betas_flat[paper_best_idx, 0],
        "best_gamma_realP": gammas_flat[realP_best_idx, 0],
        "best_beta_realP": betas_flat[realP_best_idx, 0],
        "best_gamma_realNP": gammas_flat[realNP_best_idx, 0] if realNP_best_idx is not None else None,
        "best_beta_realNP": betas_flat[realNP_best_idx, 0] if realNP_best_idx is not None else None,
        "best_gamma_qaoa": gammas_flat[real_qaoa_best_idx, 0] if real_qaoa_best_idx is not None else None,
        "best_beta_qaoa": betas_flat[real_qaoa_best_idx, 0] if real_qaoa_best_idx is not None else None,
        "paper_landscape": paper_landscape,
        "realP_landscape": realP_landscape,
        "realNP_landscape": realNP_landscape,
        "real_qaoa_landscape": real_qaoa_landscape,
        "paper_P": paper_P,
        "real_P": real_P,
        "N_real": np.asarray(N_real) if compute_real_homodist else None,
        "ising_model": ising_model,
    }


# ===== Run all seeds =====
print(f"Running experiment for {len(seeds)} graph seeds: {list(seeds)}")
print(f"  n={num_nodes}, p_edge={edge_probability}, p={p}, grid={num_gamma}x{num_beta}")
print(f"  compute_real_homodist={compute_real_homodist}, compute_real_qaoa={compute_real_qaoa}")
print(f"  compute_averaged_distributions={compute_averaged_distributions}, generate_plots={generate_plots}")

results = []
tally = {"Paper N+P": 0, "Paper N+Real P": 0}
if compute_real_homodist:
    tally["Real N+P"] = 0

for seed in seeds:
    r = run_single_seed(seed)
    results.append(r)

    # Update running tally
    candidates = {"Paper N+P": r["apx_paper"], "Paper N+Real P": r["apx_realP"]}
    if r["apx_realNP"] is not None:
        candidates["Real N+P"] = r["apx_realNP"]
    best_name = max(candidates, key=candidates.get)
    tally[best_name] += 1
    n_done = len(results)
    tally_str = ", ".join(f"{k}={v}" for k, v in tally.items())
    print(f"  Running tally ({n_done}/{len(seeds)}): {tally_str}")

# ===== Batch-averaged distributions (optional) =====
if compute_averaged_distributions:
    if not compute_real_homodist:
        print("WARNING: compute_averaged_distributions requires compute_real_homodist=True, skipping.")
        for r in results:
            r["apx_avgNP"] = None
    else:
        print(f"\n{'='*65}")
        print(f"  Computing batch-averaged N and P distributions...")
        print(f"{'='*65}")

        # Find max edges across all seeds for padding
        max_edges = max(r["num_edges"] for r in results)

        # Pad and average P distributions
        padded_Ps = np.zeros((len(results), max_edges + 1))
        for i, r in enumerate(results):
            P = r["real_P"]
            padded_Ps[i, : len(P)] = P
        avg_P = padded_Ps.mean(axis=0)

        # Pad and average N distributions: shape (num_edges+1, num_nodes+1, num_edges+1)
        padded_Ns = np.zeros((len(results), max_edges + 1, num_nodes + 1, max_edges + 1))
        for i, r in enumerate(results):
            N = r["N_real"]
            s = N.shape
            padded_Ns[i, : s[0], :, : s[2]] = N
        avg_N = padded_Ns.mean(axis=0)

        # Run proxy with averaged N
        print(f"  Running proxy with Averaged N ({num_pairs} grid points)...")
        Qs_avg = jl.QAOA_proxy_multi(avg_N, gammas_flat, betas_flat)
        Q_avg_final = Qs_avg[p]

        # Compute expectations with averaged P
        avg_expectations = jl.expectation(Q_avg_final, avg_P, num_nodes).to_numpy().flatten()
        avg_best_idx = int(np.argmax(avg_expectations))
        avg_best_gamma = gammas_flat[avg_best_idx, 0]
        avg_best_beta = betas_flat[avg_best_idx, 0]
        print(f"  Averaged proxy optimal: gamma/pi={avg_best_gamma:.4f}, beta/pi={avg_best_beta:.4f}")

        # Evaluate real QAOA at averaged-optimal params for each seed
        print(f"  Evaluating real QAOA at averaged-optimal params for each seed...")
        for r in results:
            real_at_avg = get_expectation(
                num_nodes,
                r["ising_model"],
                np.array([avg_best_gamma * np.pi]),
                np.array([avg_best_beta * np.pi]),
            )
            r["apx_avgNP"] = real_at_avg / r["c_opt"]
            r["best_gamma_avgNP"] = avg_best_gamma
            r["best_beta_avgNP"] = avg_best_beta

        apx_avgNPs_tmp = [r["apx_avgNP"] for r in results]
        print(f"  Mean Avg N+P approximation ratio: {np.mean(apx_avgNPs_tmp):.4f}")
else:
    for r in results:
        r["apx_avgNP"] = None

# ===== Summary table =====
apx_papers = [r["apx_paper"] for r in results]
apx_realPs = [r["apx_realP"] for r in results]
apx_realNPs = [r["apx_realNP"] for r in results if r["apx_realNP"] is not None]
apx_avgNPs = [r["apx_avgNP"] for r in results if r["apx_avgNP"] is not None]
apx_qaoas = [r["apx_qaoa"] for r in results if r["apx_qaoa"] is not None]
diffs_realP = [r["apx_realP"] - r["apx_paper"] for r in results]
diffs_realNP = [r["apx_realNP"] - r["apx_paper"] for r in results if r["apx_realNP"] is not None]
diffs_avgNP = [r["apx_avgNP"] - r["apx_paper"] for r in results if r["apx_avgNP"] is not None]

print(f"\n\n{'='*90}")
print(f"  SUMMARY: Effect of P(c') and N(c';d,c) on proxy parameter setting")
print(f"  G({num_nodes}, {edge_probability}), p={p}, {num_gamma}x{num_beta} grid")
print(f"{'='*90}")

# Build header and format string based on active options
hdr = f"  {'Seed':>6}  {'m':>3}  {'c_opt':>5}  {'PaperN+P':>9}  {'PaperN+RP':>9}"
sep = f"  {'-'*6}  {'-'*3}  {'-'*5}  {'-'*9}  {'-'*9}"
if compute_real_homodist:
    hdr += f"  {'RealN+P':>9}"
    sep += f"  {'-'*9}"
if compute_averaged_distributions and apx_avgNPs:
    hdr += f"  {'AvgN+P':>9}"
    sep += f"  {'-'*9}"
if compute_real_qaoa:
    hdr += f"  {'QAOA opt':>9}"
    sep += f"  {'-'*9}"
hdr += f"  {'RP-PP':>8}"
sep += f"  {'-'*8}"
if compute_real_homodist:
    hdr += f"  {'RNP-PP':>8}"
    sep += f"  {'-'*8}"
if compute_averaged_distributions and apx_avgNPs:
    hdr += f"  {'ANP-PP':>8}"
    sep += f"  {'-'*8}"
print(hdr)
print(sep)

for r in results:
    line = f"  {r['seed']:>6}  {r['num_edges']:>3}  {r['c_opt']:>5}  " f"{r['apx_paper']:>9.4f}  {r['apx_realP']:>9.4f}"
    if compute_real_homodist:
        line += f"  {r['apx_realNP']:>9.4f}"
    if compute_averaged_distributions and r.get("apx_avgNP") is not None:
        line += f"  {r['apx_avgNP']:>9.4f}"
    if compute_real_qaoa:
        line += f"  {r['apx_qaoa']:>9.4f}"
    diff_rp = r["apx_realP"] - r["apx_paper"]
    line += f"  {diff_rp:>+8.4f}"
    if compute_real_homodist:
        diff_rnp = r["apx_realNP"] - r["apx_paper"]
        line += f"  {diff_rnp:>+8.4f}"
    if compute_averaged_distributions and r.get("apx_avgNP") is not None:
        diff_anp = r["apx_avgNP"] - r["apx_paper"]
        line += f"  {diff_anp:>+8.4f}"
    print(line)

# Stats row
print(sep)
mean_line = f"  {'Mean':>6}  {'':>3}  {'':>5}  {np.mean(apx_papers):>9.4f}  {np.mean(apx_realPs):>9.4f}"
std_line = f"  {'Std':>6}  {'':>3}  {'':>5}  {np.std(apx_papers):>9.4f}  {np.std(apx_realPs):>9.4f}"
if compute_real_homodist:
    mean_line += f"  {np.mean(apx_realNPs):>9.4f}"
    std_line += f"  {np.std(apx_realNPs):>9.4f}"
if compute_averaged_distributions and apx_avgNPs:
    mean_line += f"  {np.mean(apx_avgNPs):>9.4f}"
    std_line += f"  {np.std(apx_avgNPs):>9.4f}"
if compute_real_qaoa:
    mean_line += f"  {np.mean(apx_qaoas):>9.4f}"
    std_line += f"  {np.std(apx_qaoas):>9.4f}"
mean_line += f"  {np.mean(diffs_realP):>+8.4f}"
std_line += f"  {np.std(diffs_realP):>+8.4f}"
if compute_real_homodist:
    mean_line += f"  {np.mean(diffs_realNP):>+8.4f}"
    std_line += f"  {np.std(diffs_realNP):>+8.4f}"
if compute_averaged_distributions and diffs_avgNP:
    mean_line += f"  {np.mean(diffs_avgNP):>+8.4f}"
    std_line += f"  {np.std(diffs_avgNP):>+8.4f}"
print(mean_line)
print(std_line)

# Recompute final tally including averaged distributions
final_tally = {"Paper N+P": 0, "Paper N+Real P": 0}
if compute_real_homodist:
    final_tally["Real N+P"] = 0
if compute_averaged_distributions and apx_avgNPs:
    final_tally["Avg N+P"] = 0
for r in results:
    candidates = {"Paper N+P": r["apx_paper"], "Paper N+Real P": r["apx_realP"]}
    if r.get("apx_realNP") is not None:
        candidates["Real N+P"] = r["apx_realNP"]
    if r.get("apx_avgNP") is not None:
        candidates["Avg N+P"] = r["apx_avgNP"]
    best_name = max(candidates, key=candidates.get)
    final_tally[best_name] += 1

print(f"\n  Win counts:")
tally_str = ", ".join(f"{k}: {v}/{len(seeds)}" for k, v in final_tally.items())
print(f"  {tally_str}")
print(f"{'='*90}")

# ===== Visualization (optional) =====
if generate_plots:
    import matplotlib.pyplot as plt

    num_seeds = len(seeds)
    # Columns: Paper N+P, Paper N+Real P, [Real N+P], [Real QAOA], P(c') comparison
    heatmap_cols = 2
    if compute_real_homodist:
        heatmap_cols += 1
    if compute_real_qaoa:
        heatmap_cols += 1
    num_cols = heatmap_cols + 1  # +1 for P distribution bar chart
    extent = [0, 1, 0, 0.5]

    # Figure 1: Per-seed grid (one row per seed)
    fig, axes = plt.subplots(num_seeds, num_cols, figsize=(5 * num_cols, 4 * num_seeds))
    if num_seeds == 1:
        axes = axes[np.newaxis, :]

    for row, r in enumerate(results):
        all_l = [r["paper_landscape"], r["realP_landscape"]]
        if compute_real_homodist and r["realNP_landscape"] is not None:
            all_l.append(r["realNP_landscape"])
        if compute_real_qaoa and r["real_qaoa_landscape"] is not None:
            all_l.append(r["real_qaoa_landscape"])
        vmin = min(l.min() for l in all_l)
        vmax = max(l.max() for l in all_l)

        col = 0
        heatmap_info = [
            (r["paper_landscape"], "Paper N+P", r["best_gamma_paper"], r["best_beta_paper"]),
            (r["realP_landscape"], "Paper N+Real P", r["best_gamma_realP"], r["best_beta_realP"]),
        ]
        if compute_real_homodist and r["realNP_landscape"] is not None:
            heatmap_info.append(
                (r["realNP_landscape"], "Real N+P", r["best_gamma_realNP"], r["best_beta_realNP"]),
            )
        if compute_real_qaoa and r["real_qaoa_landscape"] is not None:
            heatmap_info.append(
                (r["real_qaoa_landscape"], "Real QAOA", r["best_gamma_qaoa"], r["best_beta_qaoa"]),
            )

        for landscape, title, best_g, best_b in heatmap_info:
            ax = axes[row, col]
            im = ax.imshow(landscape, origin="lower", extent=extent, aspect="auto", cmap="viridis", vmin=vmin, vmax=vmax)
            ax.plot(best_g, best_b, "r*", markersize=12)
            if row == 0:
                ax.set_title(title, fontsize=10)
            ax.set_ylabel(f"seed={r['seed']}\nm={r['num_edges']}\n\nbeta/pi", fontsize=8)
            ax.set_xlabel("gamma/pi", fontsize=8)
            ax.tick_params(labelsize=7)
            col += 1

        # P distribution comparison (last column)
        ax_p = axes[row, num_cols - 1]
        nc = max(len(r["paper_P"]), len(r["real_P"]))
        cost_vals = np.arange(nc)
        pp = np.pad(r["paper_P"], (0, nc - len(r["paper_P"])))
        rp = np.pad(r["real_P"], (0, nc - len(r["real_P"])))
        ax_p.bar(cost_vals - 0.15, pp, width=0.3, alpha=0.8, color="tab:blue", label="Paper P")
        ax_p.bar(cost_vals + 0.15, rp, width=0.3, alpha=0.8, color="tab:orange", label="Real P")
        if row == 0:
            ax_p.set_title("P(c') Comparison", fontsize=10)
            ax_p.legend(fontsize=7)
        ax_p.set_xlabel("c'", fontsize=8)
        ax_p.tick_params(labelsize=7)

    fig.suptitle(
        f"Effect of P(c') on Proxy Landscape — G({num_nodes}, {edge_probability}), p={p}, {len(seeds)} seeds",
        fontsize=13,
    )
    plt.tight_layout()
    plt.savefig("investigate_P_distribution_all_seeds.png", dpi=120, bbox_inches="tight")
    plt.show()
    print("\nPer-seed figure saved to investigate_P_distribution_all_seeds.png")

    # Figure 2: Bar chart comparing approx ratios across seeds
    fig2, ax2 = plt.subplots(figsize=(10, 5))
    x = np.arange(num_seeds)
    has_avg = compute_averaged_distributions and bool(apx_avgNPs)
    n_bars = 2 + int(compute_real_homodist) + int(has_avg) + int(compute_real_qaoa)
    width = 0.8 / n_bars
    offset = 0
    ax2.bar(x + offset * width, apx_papers, width, label="Paper N+P", color="tab:blue", alpha=0.85)
    offset += 1
    ax2.bar(x + offset * width, apx_realPs, width, label="Paper N+Real P", color="tab:orange", alpha=0.85)
    offset += 1
    if compute_real_homodist:
        ax2.bar(x + offset * width, apx_realNPs, width, label="Real N+P", color="tab:red", alpha=0.85)
        offset += 1
    if has_avg:
        ax2.bar(x + offset * width, apx_avgNPs, width, label="Avg N+P", color="tab:purple", alpha=0.85)
        offset += 1
    if compute_real_qaoa:
        ax2.bar(x + offset * width, apx_qaoas, width, label="Real QAOA Opt", color="tab:green", alpha=0.85)
    ax2.set_xticks(x + width * (n_bars - 1) / 2)
    ax2.set_xticklabels([f"s={r['seed']}" for r in results], fontsize=7, rotation=90)
    ax2.set_ylabel("Approximation Ratio")
    mean_parts = [f"Paper N+P={np.mean(apx_papers):.4f}", f"Paper N+RP={np.mean(apx_realPs):.4f}"]
    if compute_real_homodist:
        mean_parts.append(f"Real N+P={np.mean(apx_realNPs):.4f}")
    if has_avg:
        mean_parts.append(f"Avg N+P={np.mean(apx_avgNPs):.4f}")
    if compute_real_qaoa:
        mean_parts.append(f"QAOA Opt={np.mean(apx_qaoas):.4f}")
    ax2.set_title(f"Approx Ratio Comparison — G({num_nodes}, {edge_probability}), p={p}\nMean: {', '.join(mean_parts)}")
    ax2.legend()
    ax2.set_ylim(0.5, 1.0)
    plt.tight_layout()
    plt.savefig("investigate_P_distribution_comparison.png", dpi=150, bbox_inches="tight")
    plt.show()
    print("Comparison bar chart saved to investigate_P_distribution_comparison.png")
