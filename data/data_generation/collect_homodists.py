#!/usr/bin/env python3
import grips, numpy as np, networkx as nx, argparse, os, time
from datetime import datetime
import cProfile, pstats

# NOTE: don't use ipython with this, it will break things

def make_graph(graph_type, num_nodes, seed=0, probability=0.5, ws_num_neighbors=2):
    if graph_type == 'ErdosRenyi': 
        edge_probability = probability
        graph = nx.erdos_renyi_graph(num_nodes, edge_probability, seed)
    elif graph_type == 'BarabasiAlbert':
        num_edges = max(1, int(probability * (num_nodes - 1)))  # Ensure at least one edge
        graph = nx.barabasi_albert_graph(num_nodes, num_edges, seed)
    elif graph_type == 'WattsStrogatz':
        # Watts-Strogatz graph requires parameters for number of nearest neighbors and rewiring probability
        k = max(ws_num_neighbors, int(edge_probability * (num_nodes - 1)))  # Ensure at least two neighbors
        p = probability  # Rewiring probability
        graph = nx.watts_strogatz_graph(num_nodes, k, p)
    else:
        raise ValueError(f"Unknown graph type: {graph_type}")

    return graph



def args_to_str(args, pair_seperator="_"):
    parts = []
    for k, v in sorted(vars(args).items()):  # sort for deterministic order
        parts.append(f"{k}={v}")
    return pair_seperator.join(parts)



def main(args):
    print("Running script with the following parameters:") 
    print("\t", args_to_str(args, "\n\t"))

    print("Getting graphs ...")
    seeds = range(args.seedstart, args.seedstart + args.graphs)
    make_graph_lambda = lambda seed : make_graph(
        args.graphtype, args.numnodes, seed, args.probability, args.neighbors
    )
    graphs = [make_graph_lambda(seed) for seed in seeds]
    print("Finished getting graphs.")

    start_datetime = datetime.now()

    print("Running get_homogeneous_distribution on small example to get rid of numba JIT compilation ...")
    start_time = time.perf_counter()
    dummy_graphs = [nx.erdos_renyi_graph(1, 0.5)]
    dummy_homogeneous_distribution = grips.get_homogeneous_distribution(dummy_graphs, simulator_name=args.backend)
    end_time = time.perf_counter()
    time_elapsed = end_time - start_time
    print(f"Finished in {time_elapsed:.8f} seconds")


    print("\nStarting computation at: ", start_datetime.strftime("%Y-%m-%d %H:%M:%S"))
    start_time = time.perf_counter()

    homogeneous_distribution = grips.get_homogeneous_distribution(graphs, simulator_name=args.backend)
    #with cProfile.Profile() as pr:
    #    homogeneous_distribution = grips.get_homogeneous_distribution(graphs, simulator_name=args.backend)
    #stats = pstats.Stats(pr)
    #stats.sort_stats("cumtime").print_stats(30)

    end_time = time.perf_counter()
    time_elapsed = end_time - start_time
    end_datetime = datetime.now()
    print("Finished comutation at: ", end_datetime.strftime("%Y-%m-%d %H:%M:%S"))
    print(f"Task took {time_elapsed:.8f} seconds.\n\n")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.join(script_dir, "results")
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, args_to_str(args))
    np.save(outfile, homogeneous_distribution, allow_pickle=False)

    return homogeneous_distribution



if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Collects the homogeneous distributiongfor a number of graphs of the same type, with user-specified parameters.",
    )
    parser.add_argument("graphtype", type=str, choices=["ErdosRenyi", "BarabasiAlbert", "WattsStrogatz"],
                        help="Type of graph to generate.")
    parser.add_argument("numnodes", type=int, help="Number of nodes/vertices in the graph.")
    parser.add_argument("-s", "--seedstart", default=0, type=int, help="Start of the range of seeds.")
    parser.add_argument("-g", "--graphs", default=100, type=int, help="Number of graphs to use (seeds will be contiguous range).")
    parser.add_argument("-p", "--probability", type=float, default=0.5, help="Probability of an edge between each pair of vertices (for erdos_renyi graphs) or rewiring probability (for watts_strogatz). For barabasi_albert, number of edges to attach to new nodes is `probability * (num_nodes - 1)`.")
    parser.add_argument("-n", "--neighbors", default=2, type=int, help="ws_num_neighbors (for wattz_strogatz only).")
    parser.add_argument("-b", "--backend", default="auto", choices=["auto", "python", "c", "gpu", "gpumpi"], type=str, help="Backend to use for computing maxcut costs.")

    try:
        args = parser.parse_args()
        main(args)

    except SystemExit as e:
        print("\nArgument parser failed. Did you provide the correct args?")
        parser.print_help()
