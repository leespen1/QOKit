"""Run repeatable GRIPS proxy studies across graph families.

This script turns the existing exploratory pieces into one reusable pipeline:
1. Generate a batch of graphs from a chosen family.
2. Compute the averaged homogeneous distribution over that batch.
3. Fit TriangleProxy and/or NormalProxy to the averaged distribution.
4. Optimize proxy parameters for a chosen QAOA depth and schedule type.
5. Score the resulting parameters on real QAOA across the same graph batch.

For Erdos-Renyi graphs, an optional PaperProxy calibration row can be included.
"""

from __future__ import annotations

import argparse
import os
import time
from dataclasses import dataclass

import networkx as nx
import numpy as np
import pandas as pd
import qokit.maxcut as mc
import scipy.optimize

from grips.QAOA_proxy_interface import QAOA_proxy_expectation_from_gamma_beta
from grips.QAOA_simulator import get_expectation
from grips.normal_proxy import NormalProxy
from grips.paper_proxy import PaperProxy
from grips.real_distribution import (
    distribution_mean_squared_error,
    get_homogeneous_distribution,
    get_homogeneous_distribution_from_proxy,
    normalize_homodist_slices,
)
from grips.sendai_opt import fit_proxy_to_real
from grips.solve_maxcut_exact import maxcut
from grips.triangle_proxy import TriangleProxy


SUPPORTED_GRAPH_TYPES = ("erdos_renyi", "barabasi_albert", "watts_strogatz")
SUPPORTED_PROXY_TYPES = ("triangle", "normal", "paper")
SUPPORTED_SCHEDULE_TYPES = ("full", "linear_ramp")
SUPPORTED_P_SOURCES = ("native", "empirical")


@dataclass
class ProxyFitResult:
    proxy_name: str
    proxy: object
    fitted_params: list[float] | None
    fit_mse: float | None
    pearson_mean: float | None
    pearson_min: float | None
    pearson_max: float | None


class CompositeProxy:
    """Proxy wrapper that uses N(c';d,c) from one proxy but overrides P(c') with an
    externally supplied distribution.  This lets us test the effect of swapping
    the P(c') source independently of the N source."""

    def __init__(self, n_proxy, p_distribution: np.ndarray):
        self.n_proxy = n_proxy
        self.num_constraints = n_proxy.num_constraints
        self.num_qubits = n_proxy.num_qubits
        self._p_distribution = np.asarray(p_distribution, dtype=float)

    def P_cost_distribution(self, cost: int) -> float:
        if 0 <= cost < len(self._p_distribution):
            return float(self._p_distribution[cost])
        return 0.0

    def N_cost_distribution(self, cost: int) -> float:
        return self.P_cost_distribution(cost) * (1 << self.num_qubits)

    def N_cost_distance_distribution(self, cost_1: int, distance: int, cost_2: int) -> float:
        return self.n_proxy.N_cost_distance_distribution(cost_1, distance, cost_2)


def make_graph(
    graph_type: str,
    num_nodes: int,
    seed: int,
    probability: float,
    ws_num_neighbors: int,
) -> nx.Graph:
    if graph_type == "erdos_renyi":
        graph = nx.erdos_renyi_graph(num_nodes, probability, seed=seed)
    elif graph_type == "barabasi_albert":
        num_edges = max(1, int(round(probability * (num_nodes - 1))))
        graph = nx.barabasi_albert_graph(num_nodes, num_edges, seed=seed)
    elif graph_type == "watts_strogatz":
        num_neighbors = max(ws_num_neighbors, int(round(probability * (num_nodes - 1))))
        if num_neighbors >= num_nodes:
            num_neighbors = num_nodes - 1
        if num_neighbors % 2 == 1:
            num_neighbors -= 1
        num_neighbors = max(2, num_neighbors)
        graph = nx.watts_strogatz_graph(num_nodes, num_neighbors, probability, seed=seed)
    else:
        raise ValueError(f"Unsupported graph_type: {graph_type}")

    if graph.number_of_edges() == 0 and num_nodes >= 2:
        graph.add_edge(0, 1)
    return graph


def build_graph_batch(
    graph_type: str,
    num_nodes: int,
    num_graphs: int,
    seed_start: int,
    probability: float,
    ws_num_neighbors: int,
) -> list[nx.Graph]:
    return [make_graph(graph_type, num_nodes, seed_start + offset, probability, ws_num_neighbors) for offset in range(num_graphs)]


def compute_empirical_cost_distribution(graphs: list[nx.Graph], num_constraints: int) -> np.ndarray:
    """Compute the empirical P(c') averaged over a batch of graphs.

    Enumerates all 2^n bitstrings for each graph, computes the MaxCut cost of
    each, and returns the normalised histogram (probability) of length
    num_constraints + 1.
    """
    n = graphs[0].number_of_nodes()
    num_states = 1 << n
    num_costs = num_constraints + 1
    total_P = np.zeros(num_costs, dtype=float)

    for graph in graphs:
        edges = list(graph.edges())
        if not edges:
            total_P[0] += 1.0
            continue
        edges_arr = np.array(edges, dtype=int)
        all_states = np.arange(num_states, dtype=int)
        bits = ((all_states[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int8)
        costs = np.sum(bits[:, edges_arr[:, 0]] != bits[:, edges_arr[:, 1]], axis=1)
        counts = np.bincount(costs, minlength=num_costs)[:num_costs]
        total_P += counts / num_states

    total_P /= len(graphs)
    return total_P


def schedule_bounds(schedule_type: str, depth: int) -> list[tuple[float, float]]:
    if schedule_type == "full":
        return [(0.0, np.pi)] * depth + [(0.0, np.pi / 2)] * depth
    if schedule_type == "linear_ramp":
        return [(0.0, np.pi), (0.0, np.pi), (0.0, np.pi / 2), (0.0, np.pi / 2)]
    raise ValueError(f"Unsupported schedule_type: {schedule_type}")


def default_schedule_init(schedule_type: str, depth: int) -> np.ndarray:
    if schedule_type == "full":
        gamma_init = np.linspace(0.1, 0.3, depth)
        beta_init = np.linspace(0.1, 0.2, depth)
        return np.hstack([gamma_init, beta_init])
    if schedule_type == "linear_ramp":
        return np.array([0.1, 0.3, 0.1, 0.2], dtype=float)
    raise ValueError(f"Unsupported schedule_type: {schedule_type}")


def clip_schedule_params(params: np.ndarray, schedule_type: str, depth: int) -> np.ndarray:
    bounds = schedule_bounds(schedule_type, depth)
    lower = np.array([bound[0] for bound in bounds], dtype=float)
    upper = np.array([bound[1] for bound in bounds], dtype=float)
    return np.clip(np.asarray(params, dtype=float), lower, upper)


def schedule_params_to_arrays(params: np.ndarray, schedule_type: str, depth: int) -> tuple[np.ndarray, np.ndarray]:
    params = clip_schedule_params(params, schedule_type, depth)
    if schedule_type == "full":
        gammas = params[:depth]
        betas = params[depth:]
        return gammas, betas

    gamma_start, gamma_end, beta_start, beta_end = params
    gammas = np.linspace(gamma_start, gamma_end, depth)
    betas = np.linspace(beta_start, beta_end, depth)
    return gammas, betas


def make_proxy(proxy_name: str, num_constraints: int, num_qubits: int, probability: float, params: list[float] | None = None):
    if proxy_name == "triangle":
        values = [0.0, 0.0, 1.0, 1.0] if params is None else params
        return TriangleProxy(num_constraints, num_qubits, *values)
    if proxy_name == "normal":
        values = [num_constraints / 2, 1.0, 1.0] if params is None else params
        return NormalProxy(num_constraints, num_qubits, *values)
    if proxy_name == "paper":
        return PaperProxy(num_constraints, num_qubits, probability)
    raise ValueError(f"Unsupported proxy_name: {proxy_name}")


def proxy_init_and_bounds(proxy_name: str, num_constraints: int, num_qubits: int) -> tuple[list[float], list[tuple[float, float]]]:
    if proxy_name == "triangle":
        init_params = [0.0, 0.0, 1.0, 1.0]
        bounds = [
            (0.0, num_qubits**2 / 3),
            (-10.0, 10.0),
            (0.005, 2.0),
            (0.05, 2.0),
        ]
        return init_params, bounds

    if proxy_name == "normal":
        init_params = [num_constraints / 2, 1.0, 1.0]
        bounds = [
            (0.0, float(num_constraints)),
            (0.1, max(10.0, float(num_constraints))),
            (0.1, max(10.0, float(num_constraints))),
        ]
        return init_params, bounds

    raise ValueError(f"Unsupported fitted proxy_name: {proxy_name}")


def summarize_distribution_fit(proxy, homodist: np.ndarray) -> tuple[float, float, float, float]:
    predicted = get_homogeneous_distribution_from_proxy(proxy, homodist.shape[0] - 1)
    normalized_predicted = normalize_homodist_slices(predicted)
    normalized_homodist = normalize_homodist_slices(homodist)
    mse = distribution_mean_squared_error(proxy, homodist, normalize=True)

    pearsons = []
    for cost_index in range(normalized_homodist.shape[0]):
        predicted_slice = normalized_predicted[cost_index, :, :].flatten()
        homodist_slice = normalized_homodist[cost_index, :, :].flatten()
        if np.std(predicted_slice) == 0 or np.std(homodist_slice) == 0:
            continue
        pearsons.append(np.corrcoef(predicted_slice, homodist_slice)[0, 1])
    pearsons = np.asarray(pearsons, dtype=float)
    pearsons = pearsons[~np.isnan(pearsons)]
    if pearsons.size == 0:
        return mse, np.nan, np.nan, np.nan
    return mse, float(np.mean(pearsons)), float(np.min(pearsons)), float(np.max(pearsons))


def fit_proxy(
    proxy_name: str,
    homodist: np.ndarray,
    num_constraints: int,
    num_qubits: int,
    probability: float,
    fit_max_iter: int,
    fit_fail_til_shrink: int,
    fit_fail_til_end: int,
    fit_grid_size_start: int,
) -> ProxyFitResult:
    if proxy_name == "paper":
        proxy = make_proxy(proxy_name, num_constraints, num_qubits, probability)
        mse, pearson_mean, pearson_min, pearson_max = summarize_distribution_fit(proxy, homodist)
        return ProxyFitResult(proxy_name, proxy, None, mse, pearson_mean, pearson_min, pearson_max)

    init_params, bounds = proxy_init_and_bounds(proxy_name, num_constraints, num_qubits)
    proxy = make_proxy(proxy_name, num_constraints, num_qubits, probability, init_params)
    fitted_params, fit_mse = fit_proxy_to_real(
        proxy,
        homodist,
        init_params,
        bounds,
        max_iter=fit_max_iter,
        fail_til_shrink=fit_fail_til_shrink,
        fail_til_end=fit_fail_til_end,
        grid_size_start=fit_grid_size_start,
    )
    mse, pearson_mean, pearson_min, pearson_max = summarize_distribution_fit(proxy, homodist)
    return ProxyFitResult(
        proxy_name=proxy_name,
        proxy=proxy,
        fitted_params=np.asarray(fitted_params, dtype=float).tolist(),
        fit_mse=fit_mse,
        pearson_mean=pearson_mean,
        pearson_min=pearson_min,
        pearson_max=pearson_max,
    )


def optimize_proxy_schedule(
    proxy,
    depth: int,
    schedule_type: str,
    optimizer_method: str,
    optimizer_maxiter: int,
) -> dict:
    init_params = default_schedule_init(schedule_type, depth)

    def objective(raw_params: np.ndarray) -> float:
        gammas, betas = schedule_params_to_arrays(raw_params, schedule_type, depth)
        return -QAOA_proxy_expectation_from_gamma_beta(proxy, gammas, betas)

    result = scipy.optimize.minimize(
        objective,
        init_params,
        method=optimizer_method,
        options={"maxiter": optimizer_maxiter},
    )
    best_schedule_params = clip_schedule_params(result.x, schedule_type, depth)
    gammas, betas = schedule_params_to_arrays(best_schedule_params, schedule_type, depth)
    return {
        "schedule_params": best_schedule_params,
        "gammas": gammas,
        "betas": betas,
        "proxy_expectation": -float(result.fun),
        "num_proxy_calls": int(result.nfev),
        "proxy_opt_success": bool(result.success),
        "proxy_opt_message": str(result.message),
    }


def evaluate_real_qaoa(graphs: list[nx.Graph], gammas: np.ndarray, betas: np.ndarray, simulator_name: str) -> dict:
    expectations = []
    approx_ratios = []
    optimal_values = []

    for graph in graphs:
        ising_model = mc.get_maxcut_terms(graph)
        expectation = get_expectation(graph.number_of_nodes(), ising_model, gammas, betas, simulator_name=simulator_name)
        optimum, _ = maxcut(graph)
        expectations.append(float(expectation))
        optimal_values.append(float(optimum))
        approx_ratios.append(float(expectation / optimum))

    return {
        "mean_real_expectation": float(np.mean(expectations)),
        "std_real_expectation": float(np.std(expectations)),
        "mean_approx_ratio": float(np.mean(approx_ratios)),
        "std_approx_ratio": float(np.std(approx_ratios)),
        "mean_optimum": float(np.mean(optimal_values)),
    }


def build_summary_row(
    args,
    graphs: list[nx.Graph],
    homodist: np.ndarray,
    fit_result: ProxyFitResult,
    optimization_result: dict,
    real_qaoa_metrics: dict,
) -> dict:
    mean_num_edges = float(np.mean([graph.number_of_edges() for graph in graphs]))
    max_num_edges = max(graph.number_of_edges() for graph in graphs)
    schedule_params = optimization_result["schedule_params"]

    row = {
        "graph_type": args.graph_type,
        "num_nodes": args.num_nodes,
        "num_graphs": args.num_graphs,
        "seed_start": args.seed_start,
        "edge_probability": args.edge_probability,
        "ws_num_neighbors": args.ws_num_neighbors if args.graph_type == "watts_strogatz" else None,
        "proxy_name": fit_result.proxy_name,
        "depth": args.depth,
        "schedule_type": args.schedule_type,
        "mean_num_edges": mean_num_edges,
        "max_num_edges": max_num_edges,
        "homodist_shape": "x".join(str(dim) for dim in homodist.shape),
        "fit_mse": fit_result.fit_mse,
        "pearson_mean": fit_result.pearson_mean,
        "pearson_min": fit_result.pearson_min,
        "pearson_max": fit_result.pearson_max,
        "fitted_params": None if fit_result.fitted_params is None else repr(fit_result.fitted_params),
        "schedule_params": repr(schedule_params.tolist()),
        "gammas": repr(optimization_result["gammas"].tolist()),
        "betas": repr(optimization_result["betas"].tolist()),
        "proxy_expectation": optimization_result["proxy_expectation"],
        "num_proxy_calls": optimization_result["num_proxy_calls"],
        "proxy_opt_success": optimization_result["proxy_opt_success"],
        "proxy_opt_message": optimization_result["proxy_opt_message"],
    }
    row.update(real_qaoa_metrics)
    return row


def parse_proxy_names(proxy_names_raw: str, include_paper: bool, graph_type: str) -> list[str]:
    proxy_names = [name.strip() for name in proxy_names_raw.split(",") if name.strip()]
    unsupported = sorted(set(proxy_names) - set(SUPPORTED_PROXY_TYPES))
    if unsupported:
        raise ValueError(f"Unsupported proxy names: {unsupported}")
    if include_paper and "paper" not in proxy_names:
        proxy_names.append("paper")
    return proxy_names


def default_output_path(args) -> str:
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(results_dir, exist_ok=True)
    edge_probability_slug = str(args.edge_probability).replace(".", "p")
    filename = f"proxy_study_{args.graph_type}_n{args.num_nodes}_g{args.num_graphs}" f"_ep{edge_probability_slug}" f"_p{args.depth}_{args.schedule_type}"
    if args.graph_type == "watts_strogatz":
        filename += f"_k{args.ws_num_neighbors}"
    filename += ".csv"
    return os.path.join(results_dir, filename)


def main(args) -> pd.DataFrame:
    total_start = time.perf_counter()
    proxy_names = parse_proxy_names(args.proxy_names, args.include_paper, args.graph_type)
    p_sources = [s.strip() for s in args.p_sources.split(",") if s.strip()]
    for ps in p_sources:
        if ps not in SUPPORTED_P_SOURCES:
            raise ValueError(f"Unsupported p_source: {ps}")

    graphs = build_graph_batch(
        graph_type=args.graph_type,
        num_nodes=args.num_nodes,
        num_graphs=args.num_graphs,
        seed_start=args.seed_start,
        probability=args.edge_probability,
        ws_num_neighbors=args.ws_num_neighbors,
    )

    print(f"Built {len(graphs)} graphs for graph_type={args.graph_type}.")
    print("Computing averaged homogeneous distribution...")
    homodist_start = time.perf_counter()
    homodist = get_homogeneous_distribution(graphs, simulator_name=args.cost_backend)
    homodist_runtime_sec = time.perf_counter() - homodist_start
    num_constraints = homodist.shape[0] - 1

    # Compute empirical cost distribution P(c') if requested
    empirical_P = None
    if "empirical" in p_sources:
        print("Computing empirical cost distribution P(c')...")
        empirical_P = compute_empirical_cost_distribution(graphs, num_constraints)

    # Phase 1: Fit all proxy N-distributions (done once per proxy type)
    fit_results = {}
    fit_runtimes = {}
    for proxy_name in proxy_names:
        print(f"\n--- Fitting {proxy_name} proxy N(c';d,c) ---")
        if proxy_name == "paper" and args.graph_type != "erdos_renyi":
            print("Using PaperProxy as heuristic baseline; " "analytical assumptions do not match this graph type.")
        fit_start = time.perf_counter()
        fit_results[proxy_name] = fit_proxy(
            proxy_name=proxy_name,
            homodist=homodist,
            num_constraints=num_constraints,
            num_qubits=args.num_nodes,
            probability=args.edge_probability,
            fit_max_iter=args.fit_max_iter,
            fit_fail_til_shrink=args.fit_fail_til_shrink,
            fit_fail_til_end=args.fit_fail_til_end,
            fit_grid_size_start=args.fit_grid_size_start,
        )
        fit_runtimes[proxy_name] = time.perf_counter() - fit_start

    # Phase 2: For every (N-source, P-source) combination, optimise and evaluate
    rows = []
    for proxy_name in proxy_names:
        fit_result = fit_results[proxy_name]
        for p_source in p_sources:
            config_label = f"{p_source}_P+{proxy_name}_N"
            print(f"\n=== Config: {config_label} ===")

            if p_source == "empirical":
                assert empirical_P is not None
                proxy = CompositeProxy(fit_result.proxy, empirical_P)
            else:
                proxy = fit_result.proxy

            proxy_opt_start = time.perf_counter()
            optimization_result = optimize_proxy_schedule(
                proxy,
                depth=args.depth,
                schedule_type=args.schedule_type,
                optimizer_method=args.optimizer_method,
                optimizer_maxiter=args.optimizer_maxiter,
            )
            proxy_opt_runtime_sec = time.perf_counter() - proxy_opt_start

            real_qaoa_start = time.perf_counter()
            real_qaoa_metrics = evaluate_real_qaoa(
                graphs,
                optimization_result["gammas"],
                optimization_result["betas"],
                simulator_name=args.real_backend,
            )
            real_qaoa_runtime_sec = time.perf_counter() - real_qaoa_start

            row = build_summary_row(args, graphs, homodist, fit_result, optimization_result, real_qaoa_metrics)
            row["p_source"] = p_source
            row["n_source"] = proxy_name
            row["config_label"] = config_label
            row["homodist_runtime_sec"] = homodist_runtime_sec
            row["fit_runtime_sec"] = fit_runtimes[proxy_name]
            row["proxy_opt_runtime_sec"] = proxy_opt_runtime_sec
            row["real_qaoa_runtime_sec"] = real_qaoa_runtime_sec
            row["proxy_row_runtime_sec"] = fit_runtimes[proxy_name] + proxy_opt_runtime_sec + real_qaoa_runtime_sec
            rows.append(row)

            print(f"  approx_ratio={row['mean_approx_ratio']:.4f}, " f"fit_mse={row['fit_mse'] if row['fit_mse'] is not None else float('nan'):.6f}")

    results = pd.DataFrame(rows)
    total_study_runtime_sec = time.perf_counter() - total_start
    results["study_runtime_sec"] = total_study_runtime_sec
    output_path = args.output or default_output_path(args)
    results.to_csv(output_path, index=False)
    print(f"\nSaved results to {output_path}")
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run a reusable proxy-fitting and real-QAOA evaluation study.",
    )
    parser.add_argument("graph_type", choices=SUPPORTED_GRAPH_TYPES)
    parser.add_argument("num_nodes", type=int)
    parser.add_argument("--num-graphs", type=int, default=10)
    parser.add_argument("--seed-start", type=int, default=0)
    parser.add_argument("--edge-probability", type=float, default=0.5)
    parser.add_argument("--ws-num-neighbors", type=int, default=2)
    parser.add_argument("--proxy-names", default="triangle,normal")
    parser.add_argument("--include-paper", dest="include_paper", action="store_true")
    parser.add_argument("--exclude-paper", dest="include_paper", action="store_false")
    parser.add_argument("--depth", type=int, default=1)
    parser.add_argument("--schedule-type", choices=SUPPORTED_SCHEDULE_TYPES, default="full")
    parser.add_argument(
        "--p-sources", default="native,empirical", help="Comma-separated P(c') sources: native (proxy built-in) and/or empirical (from graph instances)"
    )
    parser.add_argument("--optimizer-method", default="COBYLA")
    parser.add_argument("--optimizer-maxiter", type=int, default=200)
    parser.add_argument("--fit-max-iter", type=int, default=1000)
    parser.add_argument("--fit-fail-til-shrink", type=int, default=25)
    parser.add_argument("--fit-fail-til-end", type=int, default=50)
    parser.add_argument("--fit-grid-size-start", type=int, default=0)
    parser.add_argument("--cost-backend", default="auto", choices=["auto", "python", "c", "gpu", "gpumpi"])
    parser.add_argument("--real-backend", default="auto")
    parser.add_argument("--output")
    parser.set_defaults(include_paper=True)
    return parser


if __name__ == "__main__":
    main(build_parser().parse_args())
