"""Run multiple proxy-study configurations and summarize the results.

This script sits one level above run_proxy_study.py. It is intended for
small-to-medium exploratory sweeps where ease of use and runtime visibility are
more important than squeezing out maximum throughput.
"""

from __future__ import annotations

import argparse
import os
import time
from argparse import Namespace
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from run_proxy_study import (
    SUPPORTED_GRAPH_TYPES,
    SUPPORTED_SCHEDULE_TYPES,
    main as run_single_study,
)


DEFAULT_GRAPH_TYPES = "erdos_renyi,barabasi_albert,watts_strogatz"
DEFAULT_NODE_COUNTS = "5,7,9,11"
DEFAULT_DEPTHS = "1"
DEFAULT_EDGE_PROBABILITIES = "0.3,0.5"
DEFAULT_WS_NEIGHBOR_VALUES = "2,4"


def parse_csv_values(raw_value: str, cast) -> list:
    values = [value.strip() for value in raw_value.split(",") if value.strip()]
    return [cast(value) for value in values]


def make_sweep_name(args) -> str:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return args.sweep_name or f"proxy_sweep_{timestamp}"


def make_output_dirs(args) -> tuple[str, str]:
    results_root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "sweeps")
    os.makedirs(results_root, exist_ok=True)
    sweep_dir = os.path.join(results_root, make_sweep_name(args))
    os.makedirs(sweep_dir, exist_ok=True)
    per_study_dir = os.path.join(sweep_dir, "per_study")
    os.makedirs(per_study_dir, exist_ok=True)
    return sweep_dir, per_study_dir


def build_study_args(args, graph_type: str, num_nodes: int, depth: int, schedule_type: str, output_path: str) -> Namespace:
    return Namespace(
        graph_type=graph_type,
        num_nodes=num_nodes,
        num_graphs=args.num_graphs,
        seed_start=args.seed_start,
        edge_probability=args.edge_probability,
        ws_num_neighbors=args.ws_num_neighbors,
        proxy_names=args.proxy_names,
        include_paper=args.include_paper,
        p_sources=args.p_sources,
        depth=depth,
        schedule_type=schedule_type,
        optimizer_method=args.optimizer_method,
        optimizer_maxiter=args.optimizer_maxiter,
        fit_max_iter=args.fit_max_iter,
        fit_fail_til_shrink=args.fit_fail_til_shrink,
        fit_fail_til_end=args.fit_fail_til_end,
        fit_grid_size_start=args.fit_grid_size_start,
        cost_backend=args.cost_backend,
        real_backend=args.real_backend,
        output=output_path,
    )


def study_slug(graph_type: str, num_nodes: int, depth: int, schedule_type: str, edge_probability: float, ws_num_neighbors: int) -> str:
    edge_probability_slug = str(edge_probability).replace(".", "p")
    parts = [graph_type, f"n{num_nodes}", f"ep{edge_probability_slug}", f"p{depth}", schedule_type]
    if graph_type == "watts_strogatz":
        parts.insert(3, f"k{ws_num_neighbors}")
    return "_".join(parts)


def iter_graph_parameter_sets(args, graph_type: str) -> list[tuple[float, int]]:
    edge_probabilities = parse_csv_values(args.edge_probabilities, float)
    ws_neighbor_values = parse_csv_values(args.ws_neighbor_values, int)

    if graph_type == "watts_strogatz":
        return [(edge_probability, ws_num_neighbors) for edge_probability in edge_probabilities for ws_num_neighbors in ws_neighbor_values]

    return [(edge_probability, args.ws_num_neighbors) for edge_probability in edge_probabilities]


def parameter_tuple(row: pd.Series) -> tuple:
    if row["graph_type"] == "watts_strogatz":
        return (float(row["edge_probability"]), int(row["ws_num_neighbors"]))
    return (float(row["edge_probability"]),)


def parameter_label(graph_type: str, parameter_values: tuple) -> str:
    if graph_type == "watts_strogatz":
        edge_probability, ws_num_neighbors = parameter_values
        return f"ep={edge_probability}, k={ws_num_neighbors}"
    return f"ep={parameter_values[0]}"


def proxy_bar_label(row: pd.Series) -> str:
    p_src = row.get("p_source", "native")
    return f"{row['proxy_name']} ({p_src} P), p={int(row['depth'])}"


def make_metric_bar_figures(
    df: pd.DataFrame,
    metric: str,
    title_prefix: str,
    ylabel: str,
    output_dir: str,
    filename_prefix: str,
    error_metric: str | None = None,
) -> None:
    for graph_type in list(dict.fromkeys(df["graph_type"].tolist())):
        graph_df = df[df["graph_type"] == graph_type].copy()
        graph_df["parameter_tuple"] = graph_df.apply(parameter_tuple, axis=1)
        parameter_values = sorted(graph_df["parameter_tuple"].unique())

        num_subplots = len(parameter_values)
        num_cols = min(2, num_subplots)
        num_rows = int(np.ceil(num_subplots / num_cols))
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(7 * num_cols, 4.5 * num_rows), squeeze=False)
        axes_flat = axes.flatten()

        for axis, parameter_value in zip(axes_flat, parameter_values):
            subset = graph_df[graph_df["parameter_tuple"] == parameter_value].copy()
            subset = subset.sort_values(["num_nodes", "proxy_name", "depth", "schedule_type"])
            node_values = sorted(subset["num_nodes"].unique())

            label_rows = (
                subset[[c for c in ["proxy_name", "p_source", "depth", "schedule_type"] if c in subset.columns]]
                .drop_duplicates()
                .sort_values([c for c in ["proxy_name", "p_source", "depth", "schedule_type"] if c in subset.columns])
            )
            bar_labels = [proxy_bar_label(row) for _, row in label_rows.iterrows()]
            bar_count = len(bar_labels)
            x_positions = np.arange(len(node_values), dtype=float)
            bar_width = 0.8 / max(bar_count, 1)

            for label_index, (_, label_row) in enumerate(label_rows.iterrows()):
                mask = (
                    (subset["proxy_name"] == label_row["proxy_name"])
                    & (subset["depth"] == label_row["depth"])
                    & (subset["schedule_type"] == label_row["schedule_type"])
                )
                if "p_source" in label_row.index:
                    mask = mask & (subset["p_source"] == label_row["p_source"])
                label_subset = subset[mask]
                label_subset = label_subset.set_index("num_nodes")

                metric_values = [float(label_subset.loc[num_nodes, metric]) for num_nodes in node_values]
                if error_metric is not None:
                    error_values = [float(label_subset.loc[num_nodes, error_metric]) for num_nodes in node_values]
                else:
                    error_values = None

                offsets = x_positions - 0.4 + bar_width / 2 + label_index * bar_width
                axis.bar(
                    offsets,
                    metric_values,
                    width=bar_width,
                    label=bar_labels[label_index],
                    yerr=error_values,
                    capsize=4 if error_values is not None else 0,
                )

            axis.set_title(parameter_label(graph_type, parameter_value))
            axis.set_xlabel("num_nodes")
            axis.set_ylabel(ylabel)
            axis.set_xticks(x_positions)
            axis.set_xticklabels([str(node_value) for node_value in node_values])
            axis.grid(True, axis="y", alpha=0.3)
            axis.legend(fontsize="small")

        for axis in axes_flat[num_subplots:]:
            fig.delaxes(axis)

        fig.suptitle(f"{title_prefix}: {graph_type}")
        fig.tight_layout()
        output_path = os.path.join(output_dir, f"{filename_prefix}_{graph_type}.png")
        fig.savefig(output_path, dpi=160)
        plt.close(fig)


def write_text_report(df: pd.DataFrame, sweep_dir: str, total_runtime_sec: float, args) -> None:
    report_lines = [
        "QOKit Proxy Sweep Report",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"Sweep name: {os.path.basename(sweep_dir)}",
        f"Total runtime (sec): {total_runtime_sec:.3f}",
        f"Graph types: {args.graph_types}",
        f"Node counts: {args.node_counts}",
        f"Depths: {args.depths}",
        f"Schedule types: {args.schedule_types}",
        f"Edge probabilities: {args.edge_probabilities}",
        f"WS neighbor values: {args.ws_neighbor_values}",
        f"Proxy names: {args.proxy_names}",
        f"P sources: {args.p_sources}",
        f"Graphs per config: {args.num_graphs}",
        "",
        "Best mean approximation ratio by configuration:",
    ]

    grouped = df.sort_values("mean_approx_ratio", ascending=False).groupby(
        ["graph_type", "num_nodes", "depth", "schedule_type"], as_index=False
    )
    best_rows = grouped.head(1)
    for _, row in best_rows.iterrows():
        config = row.get("config_label", row["proxy_name"])
        report_lines.append(
            f"- {row['graph_type']}, n={row['num_nodes']}, ep={row['edge_probability']}, "
            f"k={row['ws_num_neighbors']}, p={row['depth']}, {row['schedule_type']}: "
            f"best={config} approx={row['mean_approx_ratio']:.4f} runtime={row['study_runtime_sec']:.3f}s"
        )

    report_lines.extend(
        [
            "",
            "Fastest study rows:",
        ]
    )
    fastest_rows = df.nsmallest(min(5, len(df)), "study_runtime_sec")
    for _, row in fastest_rows.iterrows():
        config = row.get("config_label", row["proxy_name"])
        report_lines.append(
            f"- {row['graph_type']}, n={row['num_nodes']}, ep={row['edge_probability']}, k={row['ws_num_neighbors']}, "
            f"config={config}, p={row['depth']}, "
            f"runtime={row['study_runtime_sec']:.3f}s, approx={row['mean_approx_ratio']:.4f}"
        )

    # P+N configuration comparison table
    if "config_label" in df.columns:
        report_lines.extend(["", "P+N Configuration Comparison:"])
        for graph_type in df["graph_type"].unique():
            for num_nodes in sorted(df[df["graph_type"] == graph_type]["num_nodes"].unique()):
                report_lines.append(f"\n  {graph_type}, n={num_nodes}:")
                sub = df[(df["graph_type"] == graph_type) & (df["num_nodes"] == num_nodes)]
                sub = sub.sort_values("mean_approx_ratio", ascending=False)
                for _, r in sub.iterrows():
                    mse_str = f"{r['fit_mse']:.6f}" if pd.notna(r.get("fit_mse")) else "N/A"
                    report_lines.append(
                        f"    {str(r['config_label']):<30s}  approx={r['mean_approx_ratio']:.4f}  "
                        f"proxy_exp={r['proxy_expectation']:.4f}  fit_mse={mse_str}"
                    )

    report_path = os.path.join(sweep_dir, "report.txt")
    with open(report_path, "w", encoding="ascii") as handle:
        handle.write("\n".join(report_lines) + "\n")


def run_sweep(args) -> pd.DataFrame:
    graph_types = parse_csv_values(args.graph_types, str)
    node_counts = parse_csv_values(args.node_counts, int)
    depths = parse_csv_values(args.depths, int)
    schedule_types = parse_csv_values(args.schedule_types, str)

    unsupported_graph_types = sorted(set(graph_types) - set(SUPPORTED_GRAPH_TYPES))
    if unsupported_graph_types:
        raise ValueError(f"Unsupported graph types: {unsupported_graph_types}")

    unsupported_schedule_types = sorted(set(schedule_types) - set(SUPPORTED_SCHEDULE_TYPES))
    if unsupported_schedule_types:
        raise ValueError(f"Unsupported schedule types: {unsupported_schedule_types}")

    if max(node_counts) > 11:
        raise ValueError("This sweep script is currently intended for 11 nodes or fewer. Increase intentionally after reviewing runtime.")

    sweep_dir, per_study_dir = make_output_dirs(args)
    summary_rows = []
    total_start = time.perf_counter()

    graph_parameter_sets = {graph_type: iter_graph_parameter_sets(args, graph_type) for graph_type in graph_types}
    total_configs = sum(len(graph_parameter_sets[graph_type]) for graph_type in graph_types) * len(node_counts) * len(depths) * len(schedule_types)
    completed = 0

    for graph_type in graph_types:
        for num_nodes in node_counts:
            for edge_probability, ws_num_neighbors in graph_parameter_sets[graph_type]:
                for depth in depths:
                    for schedule_type in schedule_types:
                        completed += 1
                        slug = study_slug(graph_type, num_nodes, depth, schedule_type, edge_probability, ws_num_neighbors)
                        output_path = os.path.join(per_study_dir, f"{slug}.csv")
                        print(
                            f"\n[{completed}/{total_configs}] Running {slug} "
                            f"with num_graphs={args.num_graphs}, proxies={args.proxy_names}"
                        )
                        study_args = build_study_args(args, graph_type, num_nodes, depth, schedule_type, output_path)
                        study_args.edge_probability = edge_probability
                        study_args.ws_num_neighbors = ws_num_neighbors
                        study_df = run_single_study(study_args)
                        study_df = study_df.copy()
                        study_df["sweep_name"] = os.path.basename(sweep_dir)
                        study_df["study_slug"] = slug
                        summary_rows.append(study_df)

    summary_df = pd.concat(summary_rows, ignore_index=True)
    summary_path = os.path.join(sweep_dir, "summary.csv")
    summary_df.to_csv(summary_path, index=False)

    total_runtime_sec = time.perf_counter() - total_start
    write_text_report(summary_df, sweep_dir, total_runtime_sec, args)

    if not args.skip_plots:
        make_metric_bar_figures(
            summary_df,
            metric="mean_approx_ratio",
            title_prefix="Mean Approximation Ratio",
            ylabel="mean_approx_ratio",
            output_dir=sweep_dir,
            filename_prefix="mean_approx_ratio_bar",
            error_metric="std_approx_ratio",
        )
        make_metric_bar_figures(
            summary_df,
            metric="fit_mse",
            title_prefix="Fit MSE",
            ylabel="fit_mse",
            output_dir=sweep_dir,
            filename_prefix="fit_mse_bar",
        )
        make_metric_bar_figures(
            summary_df,
            metric="study_runtime_sec",
            title_prefix="Study Runtime",
            ylabel="study_runtime_sec",
            output_dir=sweep_dir,
            filename_prefix="study_runtime_bar",
        )

    print(f"\nSweep complete. Summary written to {summary_path}")
    print(f"Sweep outputs stored in {sweep_dir}")
    return summary_df


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run an easy-to-use multi-configuration sweep over proxy studies.",
    )
    parser.add_argument("--graph-types", default=DEFAULT_GRAPH_TYPES)
    parser.add_argument("--node-counts", default=DEFAULT_NODE_COUNTS)
    parser.add_argument("--depths", default=DEFAULT_DEPTHS)
    parser.add_argument("--schedule-types", default="full")
    parser.add_argument("--edge-probabilities", default=DEFAULT_EDGE_PROBABILITIES)
    parser.add_argument("--ws-neighbor-values", default=DEFAULT_WS_NEIGHBOR_VALUES)
    parser.add_argument("--proxy-names", default="triangle,normal")
    parser.add_argument("--p-sources", default="native,empirical",
                        help="Comma-separated P(c') sources: native and/or empirical")
    parser.add_argument("--include-paper", "--include-paper-er", dest="include_paper", action="store_true")
    parser.add_argument("--exclude-paper", dest="include_paper", action="store_false")
    parser.add_argument("--num-graphs", type=int, default=5)
    parser.add_argument("--seed-start", type=int, default=0)
    parser.add_argument("--edge-probability", type=float, default=0.5)
    parser.add_argument("--ws-num-neighbors", type=int, default=2)
    parser.add_argument("--optimizer-method", default="COBYLA")
    parser.add_argument("--optimizer-maxiter", type=int, default=100)
    parser.add_argument("--fit-max-iter", type=int, default=250)
    parser.add_argument("--fit-fail-til-shrink", type=int, default=25)
    parser.add_argument("--fit-fail-til-end", type=int, default=50)
    parser.add_argument("--fit-grid-size-start", type=int, default=0)
    parser.add_argument("--cost-backend", default="python")
    parser.add_argument("--real-backend", default="python")
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--sweep-name")
    parser.set_defaults(include_paper=True)
    return parser


if __name__ == "__main__":
    run_sweep(build_parser().parse_args())