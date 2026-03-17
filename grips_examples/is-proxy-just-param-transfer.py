"""Study whether PaperProxy behaves like parameter transfer across ER classes.

This script sweeps over Erdős-Rényi graph *classes* and records the QAOA
parameters that are optimal according to the analytical PaperProxy. The goal is
to test whether the proxy-optimal parameters meaningfully change with the graph
class, or whether the methodology effectively collapses to one nearly fixed set
of transferable parameters.

The script supports two sweep modes:

1. ``m`` mode:
   Treat the number of constraints / edges ``m`` as the primary sweep variable.
   This is the right mode when ``n`` is fixed and you want to see how the
   PaperProxy optimum moves as the edge count changes.

2. ``prob_edge`` mode:
   Treat the Erdős-Rényi edge probability ``p_edge`` as the primary sweep
   variable. For each ``n``, the script converts the requested probability to an
   integer edge count

       m = round(p_edge * n * (n - 1) / 2)

   and builds the PaperProxy using that derived ``m``. This is the preferred
   mode for comparing across multiple values of ``n`` because it keeps the graph
   density comparable across sizes.

For each sweep case, the script:

1. Builds a ``PaperProxy`` for the selected ER class.
2. Optimizes the full QAOA schedule ``(gamma_1, ..., gamma_p, beta_1, ..., beta_p)``
   over the specified parameter ranges.
3. Stores the optimized parameters and metadata in a CSV file.
4. Optionally computes depth-1 proxy landscapes as heatmaps.
5. Produces a density / scatter style plot showing where the optimized layer
   parameters fall in the ``(gamma, beta)`` plane.

Outputs
-------

The script writes three outputs into ``OUTPUT_DIR`` or the directory passed via
``--output-dir``:

1. ``param_transfer_study_n*.csv``
   One row per sweep case. Includes:
   - ``n``: number of qubits / nodes
   - ``m``: integer edge count used to instantiate the PaperProxy
   - ``requested_prob_edge``: requested ER probability in ``prob_edge`` mode
   - ``prob_edge``: effective probability implied by ``m`` and ``n``
   - ``depth``: QAOA depth
   - optimization diagnostics and runtime
   - full parameter vectors ``gammas`` and ``betas``
   - per-layer columns ``gamma_1``, ``beta_1``, ``gamma_2``, ``beta_2``, ...

2. ``landscape_p1_n*.png``
   Heatmaps of the PaperProxy objective for depth 1. Each subplot corresponds
   to one sweep case with depth 1. The optimized point is marked with a star.

3. ``optimal_params_density_n*.png``
   Scatter / density summary of all optimized parameter pairs across the sweep.
   Points are labeled by layer index. Marker shape distinguishes ``n`` and the
   color scale tracks either requested ``p_edge`` or raw ``m``.

How to use
----------

Run from the repository root, for example:

```
/home/kerger/code/QOKit/qokitvenv/bin/python grips_examples/is-proxy-just-param-transfer.py \
  --sweep-mode prob_edge \
  --n-values 8:12 \
  --depth-values 1,2 \
  --edge-prob-values 0.2,0.35,0.5,0.65 \
  --optimizer-maxiter 200 \
  --num-gamma-grid 20 \
  --num-beta-grid 20 \
  --landscape-m-max 25 \
  --output-dir grips_examples/results/param_transfer_medium
```

/home/kerger/code/QOKit/qokitvenv/bin/python grips_examples/is-proxy-just-param-transfer.py \
  --sweep-mode prob_edge \
  --n-values 8:12 \
  --depth-values 1 \
  --edge-prob-values 0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6 \
  --optimizer-maxiter 200 \
  --num-gamma-grid 20 \
  --num-beta-grid 20 \
  --landscape-m-max 25 \
  --output-dir grips_examples/results/param_transfer_medium

For a fixed-``n`` raw-``m`` study:

```
/home/kerger/code/QOKit/qokitvenv/bin/python grips_examples/is-proxy-just-param-transfer.py \
  --sweep-mode m \
  --n-qubits 10 \
  --depth-values 1,2 \
  --m-values 5,10,15,20,25,30,35,40
```

Option reference
----------------

``--n-values``
    Grid of node counts ``n`` to sweep. Accepts either comma-separated values
    like ``8,10,12`` or inclusive colon syntax like ``8:12`` or ``8:12:2``.
    Used primarily for cross-``n`` studies.

``--n-qubits``
    Convenience alias for specifying ``n`` values directly on the command line.
    It accepts the same grid syntax as ``--n-values``. Examples:
    ``--n-qubits 10`` or ``--n-qubits 8,10,12`` or ``--n-qubits 8:12:2``.
    If provided, it takes precedence over ``--n-values``.

``--depth-values``
    QAOA depths to optimize. Accepts comma-separated or inclusive colon syntax.
    Example: ``1,2`` or ``1:3``.

``--sweep-mode``
    Chooses whether the sweep grid is driven by raw edge counts ``m`` or by ER
    probabilities ``prob_edge``. Valid values:
    - ``m``: use ``--m-values``
    - ``prob_edge``: use ``--edge-prob-values``

``--m-values``
    Grid of edge counts used in ``m`` mode. Values larger than the maximum
    possible edge count for a chosen ``n`` are automatically discarded.

``--edge-prob-values``
    Grid of ER probabilities used in ``prob_edge`` mode. Accepts
    comma-separated values like ``0.2,0.35,0.5`` or colon syntax like
    ``0.2:0.8:0.1``.

``--optimizer-method``
    SciPy optimization method used for PaperProxy parameter optimization.
    Default is ``Powell``. Other methods may work, but some behave poorly with
    the proxy objective unless bounds are well supported.

``--optimizer-maxiter``
    Maximum number of optimizer iterations. Increase this for more reliable
    optima, especially at depth 2 or above.

``--gamma-range``
    Comma-separated lower and upper bounds for gamma values, used both for the
    optimizer bounds and the plotted landscape range.

``--beta-range``
    Comma-separated lower and upper bounds for beta values, used both for the
    optimizer bounds and the plotted landscape range.

``--num-gamma-grid``
    Number of gamma grid points used when building depth-1 heatmaps.

``--num-beta-grid``
    Number of beta grid points used when building depth-1 heatmaps.

``--landscape-m-max``
    Optional cutoff for landscape generation. If provided, depth-1 heatmaps are
    only computed for cases with ``m <= landscape_m_max``. This is useful when
    large-``m`` cases are too slow to visualize densely.

``--output-dir``
    Directory where the CSV and plots will be saved.

Notes
-----

- In multi-``n`` studies, prefer ``--sweep-mode prob_edge`` rather than raw
  ``m``. Using fixed ``m`` across changing ``n`` mixes density effects with
  size effects and is harder to interpret.
- The PaperProxy is a class-level ER model. It does not encode the detailed
  topology of an individual graph instance.
- Depth-1 landscapes can be slow for large ``m`` because the PaperProxy cost
  kernel becomes expensive to evaluate repeatedly.
"""

import argparse
import ast
import itertools
import math
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import qokit.maxcut as mc
import scipy.optimize

from grips.QAOA_proxy_interface import QAOA_proxy_expectation_from_gamma_beta
from grips.QAOA_simulator import get_expectation
from grips.paper_proxy import PaperProxy
from grips.solve_maxcut_exact import maxcut


N_QUBITS = 10
N_VALUES = [N_QUBITS]
P_VALUES = [1, 2, 3]
M_VALUES = list(range(5, 46, 5))
EDGE_PROB_VALUES = [0.2, 0.35, 0.5, 0.65]
SWEEP_MODE = "m"
OPTIMIZER_METHOD = "Powell"
OPTIMIZER_MAXITER = 500
GAMMA_PLOT_RANGE = (0.0, math.pi)
BETA_PLOT_RANGE = (0.0, math.pi / 2)
NUM_GAMMA_GRID = 40
NUM_BETA_GRID = 40
LANDSCAPE_M_MAX = None
NUM_EVAL_GRAPHS = 10
EVAL_SEED_START = 1000
REAL_BACKEND = "python"
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "results")


def parse_grid_spec(spec: str) -> list[int]:
    spec = spec.strip()
    if spec and "," not in spec and ":" not in spec:
        return [int(spec)]

    if "," in spec:
        values = [int(piece.strip()) for piece in spec.split(",") if piece.strip()]
        if not values:
            raise ValueError(f"Grid spec produced no values: {spec}")
        return values

    parts = [piece.strip() for piece in spec.split(":") if piece.strip()]
    if len(parts) not in (2, 3):
        raise ValueError(f"Unsupported grid spec: {spec}")

    start = int(parts[0])
    stop = int(parts[1])
    step = int(parts[2]) if len(parts) == 3 else 1
    if step == 0:
        raise ValueError("Grid step cannot be zero.")

    if step > 0:
        stop_adjusted = stop + 1
    else:
        stop_adjusted = stop - 1
    return list(range(start, stop_adjusted, step))


def parse_float_pair(spec: str) -> tuple[float, float]:
    parts = [float(piece.strip()) for piece in spec.split(",") if piece.strip()]
    if len(parts) != 2:
        raise ValueError(f"Expected two comma-separated floats, got: {spec}")
    lower, upper = parts
    if lower >= upper:
        raise ValueError(f"Expected lower < upper in range spec, got: {spec}")
    return lower, upper


def parse_float_grid_spec(spec: str) -> list[float]:
    spec = spec.strip()
    if spec and "," not in spec and ":" not in spec:
        return [float(spec)]

    if "," in spec:
        values = [float(piece.strip()) for piece in spec.split(",") if piece.strip()]
        if not values:
            raise ValueError(f"Grid spec produced no values: {spec}")
        return values

    parts = [piece.strip() for piece in spec.split(":") if piece.strip()]
    if len(parts) not in (2, 3):
        raise ValueError(f"Unsupported grid spec: {spec}")

    start = float(parts[0])
    stop = float(parts[1])
    step = float(parts[2]) if len(parts) == 3 else 1.0
    if step == 0:
        raise ValueError("Grid step cannot be zero.")

    values = []
    current = start
    epsilon = abs(step) / 1000.0
    if step > 0:
        while current <= stop + epsilon:
            values.append(float(current))
            current += step
    else:
        while current >= stop - epsilon:
            values.append(float(current))
            current += step
    return values


def max_edges(num_qubits: int) -> int:
    return num_qubits * (num_qubits - 1) // 2


def sanitize_m_values(num_qubits: int, m_values: list[int]) -> list[int]:
    max_m = max_edges(num_qubits)
    valid_values = sorted({value for value in m_values if 0 <= value <= max_m})
    if not valid_values:
        raise ValueError(f"No valid m values for n={num_qubits}; max allowed is {max_m}.")
    return valid_values


def edge_probability_from_m(num_qubits: int, num_constraints: int) -> float:
    denominator = num_qubits * (num_qubits - 1)
    if denominator == 0:
        raise ValueError("Need at least two qubits to define an Erdos-Renyi edge probability.")
    return (2.0 * num_constraints) / denominator


def m_from_edge_probability(num_qubits: int, probability: float) -> int:
    if not 0.0 <= probability <= 1.0:
        raise ValueError(f"Edge probability must be in [0, 1], got {probability}.")
    return int(round(probability * max_edges(num_qubits)))


def sweep_slug_from_values(values: list[int | float]) -> str:
    if not values:
        return "empty"
    unique_values = list(dict.fromkeys(values))
    if len(unique_values) == 1:
        return str(unique_values[0]).replace(".", "p")
    return f"{str(unique_values[0]).replace('.', 'p')}-to-{str(unique_values[-1]).replace('.', 'p')}"


def build_sweep_cases(
    n_values: list[int],
    depth_values: list[int],
    sweep_mode: str,
    m_values: list[int] | None,
    edge_prob_values: list[float] | None,
) -> list[dict]:
    if sweep_mode not in {"m", "prob_edge"}:
        raise ValueError(f"Unsupported sweep mode: {sweep_mode}")

    cases = []
    for num_qubits in sorted(set(n_values)):
        if num_qubits < 2:
            raise ValueError("All n values must be at least 2.")

        if sweep_mode == "m":
            if not m_values:
                raise ValueError("m sweep mode requires m_values.")
            cleaned_m_values = sanitize_m_values(num_qubits, m_values)
            for num_constraints, depth in itertools.product(cleaned_m_values, depth_values):
                cases.append(
                    {
                        "n": num_qubits,
                        "m": num_constraints,
                        "depth": depth,
                        "requested_prob_edge": np.nan,
                        "prob_edge": edge_probability_from_m(num_qubits, num_constraints),
                        "sweep_mode": "m",
                    }
                )
            continue

        if not edge_prob_values:
            raise ValueError("prob_edge sweep mode requires edge_prob_values.")

        for requested_prob_edge, depth in itertools.product(edge_prob_values, depth_values):
            num_constraints = m_from_edge_probability(num_qubits, requested_prob_edge)
            cases.append(
                {
                    "n": num_qubits,
                    "m": num_constraints,
                    "depth": depth,
                    "requested_prob_edge": float(requested_prob_edge),
                    "prob_edge": edge_probability_from_m(num_qubits, num_constraints),
                    "sweep_mode": "prob_edge",
                }
            )

    return cases


def default_initial_params(depth: int) -> np.ndarray:
    gamma_init = np.linspace(0.1, 0.3, depth)
    beta_init = np.linspace(0.1, 0.2, depth)
    return np.hstack([gamma_init, beta_init]).astype(float)


def optimization_bounds(depth: int, gamma_range: tuple[float, float], beta_range: tuple[float, float]) -> list[tuple[float, float]]:
    return [gamma_range] * depth + [beta_range] * depth


def make_proxy(num_qubits: int, num_constraints: int) -> PaperProxy:
    probability = edge_probability_from_m(num_qubits, num_constraints)
    return PaperProxy(num_constraints, num_qubits, probability)


def evaluate_real_qaoa(
    num_qubits: int,
    prob_edge: float,
    gammas: np.ndarray,
    betas: np.ndarray,
    num_eval_graphs: int,
    eval_seed_start: int,
    real_backend: str,
) -> dict:
    """Generate ER graphs and evaluate real QAOA with the given parameters."""
    expectations = []
    approx_ratios = []
    optimal_values = []

    for i in range(num_eval_graphs):
        graph = nx.erdos_renyi_graph(num_qubits, prob_edge, seed=eval_seed_start + i)
        if graph.number_of_edges() == 0 and num_qubits >= 2:
            graph.add_edge(0, 1)
        ising_model = mc.get_maxcut_terms(graph)
        expectation = get_expectation(num_qubits, ising_model, gammas, betas, simulator_name=real_backend)
        optimum, _ = maxcut(graph)
        expectations.append(float(expectation))
        optimal_values.append(float(optimum))
        if optimum > 0:
            approx_ratios.append(float(expectation / optimum))

    return {
        "mean_real_expectation": float(np.mean(expectations)),
        "std_real_expectation": float(np.std(expectations)),
        "mean_approx_ratio": float(np.mean(approx_ratios)) if approx_ratios else float("nan"),
        "std_approx_ratio": float(np.std(approx_ratios)) if approx_ratios else float("nan"),
        "mean_optimum": float(np.mean(optimal_values)),
    }


def compute_mean_params(rows: list[dict], depth_values: list[int]) -> dict[int, dict]:
    """Compute the mean gammas and betas across all sweep cases, grouped by (n, depth)."""
    mean_params: dict[int, dict] = {}
    for depth in depth_values:
        depth_rows = [r for r in rows if r["depth"] == depth]
        if not depth_rows:
            continue
        all_gammas = np.array([r["gammas"] for r in depth_rows], dtype=float)
        all_betas = np.array([r["betas"] for r in depth_rows], dtype=float)
        mean_params[depth] = {
            "gammas": np.mean(all_gammas, axis=0),
            "betas": np.mean(all_betas, axis=0),
            "std_gammas": np.std(all_gammas, axis=0),
            "std_betas": np.std(all_betas, axis=0),
            "num_cases": len(depth_rows),
        }
    return mean_params


def optimize_proxy_parameters(
    proxy: PaperProxy,
    depth: int,
    optimizer_method: str,
    optimizer_maxiter: int,
    gamma_range: tuple[float, float],
    beta_range: tuple[float, float],
) -> dict:
    initial_params = default_initial_params(depth)
    bounds = optimization_bounds(depth, gamma_range, beta_range)

    def objective(raw_params: np.ndarray) -> float:
        gammas = np.asarray(raw_params[:depth], dtype=float)
        betas = np.asarray(raw_params[depth:], dtype=float)
        return -QAOA_proxy_expectation_from_gamma_beta(proxy, gammas, betas)

    start_time = time.perf_counter()
    result = scipy.optimize.minimize(
        objective,
        initial_params,
        method=optimizer_method,
        bounds=bounds,
        options={"maxiter": optimizer_maxiter},
    )
    runtime = time.perf_counter() - start_time

    optimal_params = np.asarray(result.x, dtype=float)
    gammas = optimal_params[:depth]
    betas = optimal_params[depth:]
    return {
        "gammas": gammas,
        "betas": betas,
        "proxy_expectation": -float(result.fun),
        "success": bool(result.success),
        "message": str(result.message),
        "num_calls": int(result.nfev),
        "runtime_sec": runtime,
    }


def build_landscape(
    proxy: PaperProxy,
    gamma_range: tuple[float, float],
    beta_range: tuple[float, float],
    num_gamma_grid: int,
    num_beta_grid: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    gamma_values = np.linspace(gamma_range[0], gamma_range[1], num_gamma_grid)
    beta_values = np.linspace(beta_range[0], beta_range[1], num_beta_grid)
    landscape = np.zeros((num_beta_grid, num_gamma_grid), dtype=float)

    for beta_index, beta in enumerate(beta_values):
        for gamma_index, gamma in enumerate(gamma_values):
            landscape[beta_index, gamma_index] = QAOA_proxy_expectation_from_gamma_beta(
                proxy,
                np.array([gamma], dtype=float),
                np.array([beta], dtype=float),
            )

    return gamma_values, beta_values, landscape


def results_to_dataframe(rows: list[dict], max_depth: int) -> pd.DataFrame:
    normalized_rows = []
    for row in rows:
        normalized = dict(row)
        for layer_index in range(max_depth):
            gamma_key = f"gamma_{layer_index + 1}"
            beta_key = f"beta_{layer_index + 1}"
            normalized[gamma_key] = row["gammas"][layer_index] if layer_index < len(row["gammas"]) else np.nan
            normalized[beta_key] = row["betas"][layer_index] if layer_index < len(row["betas"]) else np.nan
        normalized["gammas"] = repr([float(value) for value in row["gammas"]])
        normalized["betas"] = repr([float(value) for value in row["betas"]])
        normalized_rows.append(normalized)

    dataframe = pd.DataFrame(normalized_rows)
    sort_columns = [column for column in ["n", "depth", "requested_prob_edge", "m"] if column in dataframe.columns]
    return dataframe.sort_values(sort_columns).reset_index(drop=True)


def save_landscape_plot(
    landscape_rows: list[dict],
    output_path: str,
    gamma_range: tuple[float, float],
    beta_range: tuple[float, float],
) -> None:
    if not landscape_rows:
        return

    num_plots = len(landscape_rows)
    num_cols = min(3, num_plots)
    num_rows = math.ceil(num_plots / num_cols)
    fig, axes = plt.subplots(
        num_rows,
        num_cols,
        figsize=(5.5 * num_cols, 4.5 * num_rows),
        squeeze=False,
        constrained_layout=True,
    )

    active_axes = []
    image = None
    for axis, row in zip(axes.flat, landscape_rows):
        active_axes.append(axis)
        image = axis.imshow(
            row["landscape"],
            origin="lower",
            aspect="auto",
            extent=(gamma_range[0], gamma_range[1], beta_range[0], beta_range[1]),
            cmap="viridis",
        )
        axis.scatter(row["gamma_star"], row["beta_star"], marker="*", s=120, c="white", edgecolors="black", linewidths=0.8)
        title = f"n={row['n']}, m={row['m']}, p_edge={row['prob_edge']:.3f}"
        if not np.isnan(row["requested_prob_edge"]):
            title += f"\nreq={row['requested_prob_edge']:.3f}"
        axis.set_title(title)
        axis.set_xlabel("gamma")
        axis.set_ylabel("beta")

    for axis in axes.flat[num_plots:]:
        axis.remove()

    if image is not None and active_axes:
        fig.colorbar(image, ax=active_axes, label="Proxy expectation")
    fig.suptitle("PaperProxy landscapes at depth 1")
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def save_density_plot(
    dataframe: pd.DataFrame,
    depth_values: list[int],
    output_path: str,
    gamma_range: tuple[float, float],
    beta_range: tuple[float, float],
) -> None:
    if dataframe.empty:
        return

    fig, axes = plt.subplots(1, len(depth_values), figsize=(6 * len(depth_values), 5), squeeze=False, constrained_layout=True)

    for axis, depth in zip(axes.flat, depth_values):
        subset = dataframe[dataframe["depth"] == depth]
        if subset.empty:
            axis.set_xlim(*gamma_range)
            axis.set_ylim(*beta_range)
            axis.set_title(f"Optimized parameters, depth={depth}")
            axis.set_xlabel("gamma")
            axis.set_ylabel("beta")
            continue

        point_rows = []
        for _, row in subset.iterrows():
            gammas = np.array(ast.literal_eval(row["gammas"]), dtype=float)
            betas = np.array(ast.literal_eval(row["betas"]), dtype=float)
            color_value = row["requested_prob_edge"] if not np.isnan(row["requested_prob_edge"]) else row["m"]
            for layer_index, (gamma, beta) in enumerate(zip(gammas, betas), start=1):
                point_rows.append(
                    {
                        "gamma": gamma,
                        "beta": beta,
                        "n": int(row["n"]),
                        "m": int(row["m"]),
                        "layer": layer_index,
                        "color_value": float(color_value),
                    }
                )

        if point_rows:
            points = pd.DataFrame(point_rows)
            gamma_points = points["gamma"].to_numpy(dtype=float)
            beta_points = points["beta"].to_numpy(dtype=float)
            axis.hexbin(
                gamma_points,
                beta_points,
                gridsize=20,
                extent=(gamma_range[0], gamma_range[1], beta_range[0], beta_range[1]),
                mincnt=1,
                cmap="Blues",
            )

            marker_cycle = ["o", "s", "^", "D", "P", "X", "v", "<", ">"]
            scatter = None
            for marker_index, num_qubits in enumerate(sorted(points["n"].unique())):
                marker = marker_cycle[marker_index % len(marker_cycle)]
                n_subset = points[points["n"] == num_qubits]
                scatter = axis.scatter(
                    n_subset["gamma"],
                    n_subset["beta"],
                    c=n_subset["color_value"],
                    cmap="plasma",
                    s=60,
                    marker=marker,
                    edgecolors="black",
                    linewidths=0.5,
                    alpha=0.9,
                    label=f"n={num_qubits}",
                )
                for _, point in n_subset.iterrows():
                    axis.text(point["gamma"], point["beta"], str(int(point["layer"])), fontsize=7, ha="center", va="center", color="white")

            if scatter is not None:
                color_label = "requested p_edge" if subset["requested_prob_edge"].notna().any() else "m"
                fig.colorbar(scatter, ax=axis, label=color_label)
            axis.legend(loc="best", fontsize=8)

        axis.set_xlim(*gamma_range)
        axis.set_ylim(*beta_range)
        axis.set_title(f"Optimized parameters, depth={depth}")
        axis.set_xlabel("gamma")
        axis.set_ylabel("beta")

    fig.suptitle("PaperProxy optimal parameter locations")
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def save_transfer_comparison_plot(
    comparison_df: pd.DataFrame,
    depth_values: list[int],
    output_path: str,
) -> None:
    """Bar chart comparing per-case proxy params vs mean transferred params."""
    if comparison_df.empty:
        return

    for depth in depth_values:
        sub = comparison_df[comparison_df["depth"] == depth]
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=(max(8, len(sub) * 0.8), 5), constrained_layout=True)
        x = np.arange(len(sub))
        width = 0.35

        individual = sub["individual_approx_ratio"].to_numpy(dtype=float)
        transferred = sub["transferred_approx_ratio"].to_numpy(dtype=float)

        bars1 = ax.bar(x - width / 2, individual, width, label="Per-case proxy params",
                       yerr=sub["individual_std_approx_ratio"].to_numpy(dtype=float), capsize=3)
        bars2 = ax.bar(x + width / 2, transferred, width, label="Mean transferred params",
                       yerr=sub["transferred_std_approx_ratio"].to_numpy(dtype=float), capsize=3)

        labels = []
        for _, row in sub.iterrows():
            pe = row.get("requested_prob_edge")
            if pd.notna(pe):
                labels.append(f"n={int(row['n'])}\nm={int(row['m'])}\npe={pe:.2f}")
            else:
                labels.append(f"n={int(row['n'])}\nm={int(row['m'])}")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_ylabel("Mean Approx Ratio")
        ax.set_title(f"Per-case Proxy Params vs Mean Transferred Params (depth={depth})")
        ax.legend()
        ax.grid(True, axis="y", alpha=0.3)

        depth_path = output_path.replace(".png", f"_d{depth}.png")
        fig.savefig(depth_path, dpi=160, bbox_inches="tight")
        plt.close(fig)


def run_study(
    n_values: list[int],
    depth_values: list[int],
    sweep_mode: str,
    m_values: list[int] | None,
    edge_prob_values: list[float] | None,
    optimizer_method: str,
    optimizer_maxiter: int,
    gamma_range: tuple[float, float],
    beta_range: tuple[float, float],
    num_gamma_grid: int,
    num_beta_grid: int,
    landscape_m_max: int | None,
    output_dir: str,
    num_eval_graphs: int = NUM_EVAL_GRAPHS,
    eval_seed_start: int = EVAL_SEED_START,
    real_backend: str = REAL_BACKEND,
) -> dict:
    os.makedirs(output_dir, exist_ok=True)
    depth_values = sorted(set(depth_values))
    if any(depth < 1 for depth in depth_values):
        raise ValueError("All depth values must be positive integers.")
    sweep_cases = build_sweep_cases(n_values, depth_values, sweep_mode, m_values, edge_prob_values)

    rows = []
    landscape_rows = []

    for case in sweep_cases:
        num_qubits = case["n"]
        num_constraints = case["m"]
        depth = case["depth"]
        probability = case["prob_edge"]
        proxy = make_proxy(num_qubits, num_constraints)
        print(f"Optimizing PaperProxy for n={num_qubits}, m={num_constraints}, depth={depth}, p_edge={probability:.3f}")
        optimization = optimize_proxy_parameters(
            proxy,
            depth,
            optimizer_method,
            optimizer_maxiter,
            gamma_range,
            beta_range,
        )

        row = {
            "sweep_mode": case["sweep_mode"],
            "n": num_qubits,
            "m": num_constraints,
            "requested_prob_edge": case["requested_prob_edge"],
            "depth": depth,
            "prob_edge": probability,
            "proxy_expectation": optimization["proxy_expectation"],
            "classical_opt_success": optimization["success"],
            "scipy_opt_message": optimization["message"],
            "num_proxy_calls": optimization["num_calls"],
            "runtime_sec": optimization["runtime_sec"],
            "gammas": optimization["gammas"],
            "betas": optimization["betas"],
        }
        rows.append(row)

        should_plot_landscape = depth == 1 and (landscape_m_max is None or num_constraints <= landscape_m_max)
        if should_plot_landscape:
            gamma_values, beta_values, landscape = build_landscape(
                proxy,
                gamma_range,
                beta_range,
                num_gamma_grid,
                num_beta_grid,
            )
            landscape_rows.append(
                {
                    "n": num_qubits,
                    "m": num_constraints,
                    "requested_prob_edge": case["requested_prob_edge"],
                    "prob_edge": probability,
                    "gamma_values": gamma_values,
                    "beta_values": beta_values,
                    "landscape": landscape,
                    "gamma_star": optimization["gammas"][0],
                    "beta_star": optimization["betas"][0],
                }
            )

    dataframe = results_to_dataframe(rows, max(depth_values))
    n_slug = sweep_slug_from_values(sorted(set(int(case["n"]) for case in sweep_cases)))
    csv_path = os.path.join(output_dir, f"param_transfer_study_n{n_slug}.csv")
    dataframe.to_csv(csv_path, index=False)

    landscape_path = os.path.join(output_dir, f"landscape_p1_n{n_slug}.png")
    save_landscape_plot(landscape_rows, landscape_path, gamma_range, beta_range)

    density_path = os.path.join(output_dir, f"optimal_params_density_n{n_slug}.png")
    save_density_plot(dataframe, depth_values, density_path, gamma_range, beta_range)

    # --- Real QAOA evaluation: per-case params vs mean transferred params ---
    print("\n=== Real QAOA transfer comparison ===")
    comparison_rows = []

    # Group rows by n to compute mean params per (n, depth)
    n_groups: dict[int, list[dict]] = {}
    for row in rows:
        n_groups.setdefault(row["n"], []).append(row)

    for num_qubits, n_rows in sorted(n_groups.items()):
        mean_params = compute_mean_params(n_rows, depth_values)

        for row in n_rows:
            depth = row["depth"]
            prob_edge = row["prob_edge"]

            # Evaluate with this case's individual proxy-optimal params
            print(f"  Evaluating n={num_qubits}, m={row['m']}, depth={depth}, p_edge={prob_edge:.3f} -- individual params")
            individual_metrics = evaluate_real_qaoa(
                num_qubits, prob_edge, row["gammas"], row["betas"],
                num_eval_graphs, eval_seed_start, real_backend,
            )

            # Evaluate with the mean transferred params for this (n, depth)
            mp = mean_params[depth]
            print(f"  Evaluating n={num_qubits}, m={row['m']}, depth={depth}, p_edge={prob_edge:.3f} -- mean transferred params")
            transferred_metrics = evaluate_real_qaoa(
                num_qubits, prob_edge, mp["gammas"], mp["betas"],
                num_eval_graphs, eval_seed_start, real_backend,
            )

            comparison_rows.append({
                "n": num_qubits,
                "m": row["m"],
                "depth": depth,
                "requested_prob_edge": row["requested_prob_edge"],
                "prob_edge": prob_edge,
                "individual_approx_ratio": individual_metrics["mean_approx_ratio"],
                "individual_std_approx_ratio": individual_metrics["std_approx_ratio"],
                "individual_expectation": individual_metrics["mean_real_expectation"],
                "transferred_approx_ratio": transferred_metrics["mean_approx_ratio"],
                "transferred_std_approx_ratio": transferred_metrics["std_approx_ratio"],
                "transferred_expectation": transferred_metrics["mean_real_expectation"],
                "mean_optimum": individual_metrics["mean_optimum"],
                "individual_gammas": repr([float(v) for v in row["gammas"]]),
                "individual_betas": repr([float(v) for v in row["betas"]]),
                "transferred_gammas": repr([float(v) for v in mp["gammas"]]),
                "transferred_betas": repr([float(v) for v in mp["betas"]]),
                "num_cases_in_mean": mp["num_cases"],
                "num_eval_graphs": num_eval_graphs,
            })

            diff = individual_metrics["mean_approx_ratio"] - transferred_metrics["mean_approx_ratio"]
            print(f"    individual={individual_metrics['mean_approx_ratio']:.4f}, "
                  f"transferred={transferred_metrics['mean_approx_ratio']:.4f}, "
                  f"diff={diff:+.4f}")

    comparison_df = pd.DataFrame(comparison_rows)
    comparison_csv_path = os.path.join(output_dir, f"transfer_comparison_n{n_slug}.csv")
    comparison_df.to_csv(comparison_csv_path, index=False)
    print(f"\nSaved transfer comparison to {comparison_csv_path}")

    # Summary statistics
    if not comparison_df.empty:
        mean_ind = comparison_df["individual_approx_ratio"].mean()
        mean_xfer = comparison_df["transferred_approx_ratio"].mean()
        ind_wins = (comparison_df["individual_approx_ratio"] > comparison_df["transferred_approx_ratio"]).sum()
        total = len(comparison_df)
        print(f"\nTransfer comparison summary:")
        print(f"  Mean individual approx ratio: {mean_ind:.4f}")
        print(f"  Mean transferred approx ratio: {mean_xfer:.4f}")
        print(f"  Individual wins: {ind_wins}/{total}")

    transfer_plot_path = os.path.join(output_dir, f"transfer_comparison_n{n_slug}.png")
    save_transfer_comparison_plot(comparison_df, depth_values, transfer_plot_path)

    return {
        "dataframe": dataframe,
        "comparison_dataframe": comparison_df,
        "csv_path": csv_path,
        "comparison_csv_path": comparison_csv_path,
        "landscape_path": landscape_path,
        "density_path": density_path,
        "transfer_plot_path": transfer_plot_path,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Study whether PaperProxy behaves like parameter transfer across ER graph classes.")
    parser.add_argument("--n-values", type=str, default=",".join(str(value) for value in N_VALUES))
    parser.add_argument("--n-qubits", type=str, default=None)
    parser.add_argument("--depth-values", type=str, default=",".join(str(value) for value in P_VALUES))
    parser.add_argument("--sweep-mode", choices=["m", "prob_edge"], default=SWEEP_MODE)
    parser.add_argument("--m-values", type=str, default=",".join(str(value) for value in M_VALUES))
    parser.add_argument("--edge-prob-values", type=str, default=",".join(str(value) for value in EDGE_PROB_VALUES))
    parser.add_argument("--optimizer-method", type=str, default=OPTIMIZER_METHOD)
    parser.add_argument("--optimizer-maxiter", type=int, default=OPTIMIZER_MAXITER)
    parser.add_argument("--gamma-range", type=str, default=f"{GAMMA_PLOT_RANGE[0]},{GAMMA_PLOT_RANGE[1]}")
    parser.add_argument("--beta-range", type=str, default=f"{BETA_PLOT_RANGE[0]},{BETA_PLOT_RANGE[1]}")
    parser.add_argument("--num-gamma-grid", type=int, default=NUM_GAMMA_GRID)
    parser.add_argument("--num-beta-grid", type=int, default=NUM_BETA_GRID)
    parser.add_argument("--landscape-m-max", type=int, default=LANDSCAPE_M_MAX)
    parser.add_argument("--num-eval-graphs", type=int, default=NUM_EVAL_GRAPHS,
                        help="Number of random ER graphs to generate for real QAOA evaluation")
    parser.add_argument("--eval-seed-start", type=int, default=EVAL_SEED_START,
                        help="Starting seed for evaluation graph generation")
    parser.add_argument("--real-backend", type=str, default=REAL_BACKEND,
                        help="Simulator backend for real QAOA evaluation")
    parser.add_argument("--output-dir", type=str, default=OUTPUT_DIR)
    return parser


def main(argv: list[str] | None = None) -> pd.DataFrame:
    if argv is None:
        result = run_study(
            n_values=N_VALUES,
            depth_values=P_VALUES,
            sweep_mode=SWEEP_MODE,
            m_values=M_VALUES,
            edge_prob_values=EDGE_PROB_VALUES,
            optimizer_method=OPTIMIZER_METHOD,
            optimizer_maxiter=OPTIMIZER_MAXITER,
            gamma_range=GAMMA_PLOT_RANGE,
            beta_range=BETA_PLOT_RANGE,
            num_gamma_grid=NUM_GAMMA_GRID,
            num_beta_grid=NUM_BETA_GRID,
            landscape_m_max=LANDSCAPE_M_MAX,
            output_dir=OUTPUT_DIR,
        )
        print(f"Saved CSV to {result['csv_path']}")
        print(f"Saved transfer comparison to {result['comparison_csv_path']}")
        print(f"Saved landscape plot to {result['landscape_path']}")
        print(f"Saved density plot to {result['density_path']}")
        return result["dataframe"]

    parser = build_arg_parser()
    args = parser.parse_args(argv)
    n_values = parse_grid_spec(args.n_qubits) if args.n_qubits is not None else parse_grid_spec(args.n_values)
    result = run_study(
        n_values=n_values,
        depth_values=parse_grid_spec(args.depth_values),
        sweep_mode=args.sweep_mode,
        m_values=parse_grid_spec(args.m_values),
        edge_prob_values=parse_float_grid_spec(args.edge_prob_values),
        optimizer_method=args.optimizer_method,
        optimizer_maxiter=args.optimizer_maxiter,
        gamma_range=parse_float_pair(args.gamma_range),
        beta_range=parse_float_pair(args.beta_range),
        num_gamma_grid=args.num_gamma_grid,
        num_beta_grid=args.num_beta_grid,
        landscape_m_max=args.landscape_m_max,
        output_dir=args.output_dir,
        num_eval_graphs=args.num_eval_graphs,
        eval_seed_start=args.eval_seed_start,
        real_backend=args.real_backend,
    )
    print(f"Saved CSV to {result['csv_path']}")
    print(f"Saved transfer comparison to {result['comparison_csv_path']}")
    print(f"Saved landscape plot to {result['landscape_path']}")
    print(f"Saved density plot to {result['density_path']}")
    return result["dataframe"]


if __name__ == "__main__":
    main(sys.argv[1:])