import grips, argparse, networkx as nx, numpy as np
from collections import OrderedDict
from pathlib import Path
import h5py


def write_homodist_hdf5(filepath, N_dist_mean, args_dict, seeds, num_graphs):
    """Write averaged N(c'; d, c) homogeneous distribution to HDF5 file.

    Structure:
    - /N_distribution: 3D array (num_costs, num_distances, num_costs)
    - /seeds: 1D array of graph seeds used
    - Attributes: all parameters from args_dict plus metadata
    """
    with h5py.File(filepath, "w") as f:
        # Store distribution with compression
        f.create_dataset("N_distribution", data=N_dist_mean, compression="gzip", compression_opts=4)
        f.create_dataset("seeds", data=np.array(seeds))

        # Store parameters as attributes
        for key, val in args_dict.items():
            f.attrs[key] = val

        # Add metadata about the distribution
        f.attrs["distribution_format"] = "N[c', d, c] - averaged over graphs"
        f.attrs["num_graphs_averaged"] = num_graphs
        f.attrs["index_order"] = "Python (0-indexed): N_distribution[c_prime, distance, cost]"


def write_costs_text(filepath, costs_list, args_dict):
    """Write raw costs (cost of each bitstring) to human-readable text file.

    Format:
    - Header with parameters
    - One line per graph: space-separated costs for bitstrings 0,1,2,...,2^n-1
    """
    with open(filepath, "w") as f:
        # Header
        f.write(f"# Raw costs for {args_dict['graphType']} graphs\n")
        f.write(f"# Parameters: {grips.args_to_str(args_dict, ', ')}\n")
        f.write(f"# Each row: cost(x) for x=0,1,2,...,2^n-1 (space-separated)\n")
        f.write(f"# Number of graphs: {len(costs_list)}\n")
        f.write("#\n")

        for costs in costs_list:
            f.write(" ".join(map(str, costs)) + "\n")


def write_costs_hdf5(filepath, costs_list, args_dict, seeds):
    """Write raw costs (cost of each bitstring) to HDF5 file with compression.

    Structure:
    - /costs: 2D array (num_graphs x 2^n)
    - /seeds: 1D array of graph seeds
    - Attributes: all parameters from args_dict
    """
    costs_array = np.array(costs_list)

    with h5py.File(filepath, "w") as f:
        # Store data with gzip compression
        f.create_dataset("costs", data=costs_array, compression="gzip", compression_opts=4)
        f.create_dataset("seeds", data=np.array(seeds))

        # Store all parameters as attributes
        for key, val in args_dict.items():
            f.attrs[key] = val


def write_text(filepath, Nc_list, args_dict):
    """Write N(c) distributions to human-readable text file.

    Format:
    - Header with parameters
    - One line per graph: space-separated counts for c=0,1,2,...
    """
    with open(filepath, "w") as f:
        # Header
        f.write(f"# N(c) distributions for {args_dict['graphType']} graphs\n")
        f.write(f"# Parameters: {grips.args_to_str(args_dict, ', ')}\n")
        f.write(f"# Each row: N(c) for c=0,1,2,...,max_cost (space-separated)\n")
        f.write(f"# Number of graphs: {len(Nc_list)}\n")
        f.write("#\n")

        for N_c in Nc_list:
            f.write(" ".join(map(str, N_c)) + "\n")


def write_hdf5(filepath, Nc_list, args_dict, seeds):
    """Write N(c) distributions to HDF5 file with compression.

    Structure:
    - /Nc: 2D array (num_graphs x max_cost+1), padded with zeros
    - /seeds: 1D array of graph seeds
    - Attributes: all parameters from args_dict
    - Key attributes: graphType, numNodes, nearestNeighbors, rewiringProbability
    """

    # Pad arrays to same length for efficient storage
    max_len = max(len(N_c) for N_c in Nc_list)
    Nc_padded = np.zeros((len(Nc_list), max_len), dtype=np.int64)
    for i, N_c in enumerate(Nc_list):
        Nc_padded[i, :len(N_c)] = N_c

    with h5py.File(filepath, "w") as f:
        # Store data with gzip compression (good balance of speed/size)
        f.create_dataset("Nc", data=Nc_padded, compression="gzip", compression_opts=4)
        f.create_dataset("seeds", data=np.array(seeds))

        # Store important graph parameters explicitly
        f.attrs["graphType"] = args_dict["graphType"]
        f.attrs["numNodes"] = args_dict["numNodes"]
        f.attrs["nearestNeighbors"] = args_dict["nearestNeighbors"]
        f.attrs["rewiringProbability"] = args_dict["rewiringProbability"]

        # Store all other parameters as attributes
        for key, val in args_dict.items():
            if key not in f.attrs:
                f.attrs[key] = val


def main(args):
    args_dict = OrderedDict(sorted(vars(args).items()))
    args_dict.pop("backend")  # Don't keep track of backend
    args_dict.pop("format")   # Don't keep track of format
    store_costs = args_dict.pop("storeCosts")  # Don't keep track of storeCosts
    compute_homodist = args_dict.pop("computeHomodist")  # Don't keep track in filename
    graphType = "WattsStrogatz"
    args_dict["graphType"] = graphType

    # Rearrange arg order for filename, make easier to parse
    args_dict.move_to_end("rewiringProbability", last=False)
    args_dict.move_to_end("nearestNeighbors", last=False)
    args_dict.move_to_end("numNodes", last=False)
    args_dict.move_to_end("graphType", last=False)
    args_dict.move_to_end("seedStart", last=True)
    args_dict.move_to_end("numGraphs", last=True)

    print("Running script with the following parameters:")
    print("\t", grips.args_to_str(args_dict, "\n\t"))

    n = args.numNodes
    k = args.nearestNeighbors
    p = args.rewiringProbability
    seeds = list(range(args.seedStart, args.seedStart + args.numGraphs))

    # Directory of this script
    script_dir = Path(__file__).resolve().parent

    # Data directory relative to script
    data_dir = script_dir / f"Data_graphType={graphType}/numNodes={n}"
    data_dir.mkdir(exist_ok=True, parents=True)

    # Initialize Julia distribution accumulator if needed
    homodist_accumulator = None
    if compute_homodist:
        from grips.julia_distributions import HomogeneousDistributionAccumulator
        # Pre-allocate for expected max edges (n*(n-1)/2 for complete graph)
        max_possible_edges = n * (n - 1) // 2
        homodist_accumulator = HomogeneousDistributionAccumulator(max_num_edges=max_possible_edges)
        use_gpu = args.backend in ("gpu", "gpumpi")
        print(f"Computing N(c',d,c) distributions using {'GPU' if use_gpu else 'CPU'}...")

    # Compute all N(c) distributions (and optionally store raw costs / compute homodist)
    Nc_list = []
    costs_list = []
    for i, seed in enumerate(seeds):
        graph = nx.watts_strogatz_graph(n, k, p, seed=seed)
        costs = np.rint(grips.get_costs(graph, args.backend)).astype(int)
        N_c = np.bincount(costs)
        Nc_list.append(N_c)

        if store_costs:
            costs_list.append(costs)

        if compute_homodist:
            num_edges = graph.number_of_edges()
            homodist_accumulator.add(costs, num_edges, n, use_gpu=use_gpu)
            if (i + 1) % 10 == 0 or (i + 1) == len(seeds):
                print(f"  Processed {i + 1}/{len(seeds)} graphs")

    # Output file base name (no extension)
    basename = grips.args_to_str(args_dict, "_")

    # Write requested format(s)
    fmt = args.format
    if fmt in ("text", "both"):
        filepath = data_dir / (basename + ".txt")
        write_text(filepath, Nc_list, args_dict)
        print(f"Wrote text file: {filepath}")

    if fmt in ("hdf5", "both"):
        filepath = data_dir / (basename + ".h5")
        write_hdf5(filepath, Nc_list, args_dict, seeds)
        print(f"Wrote HDF5 file: {filepath}")

    # Write raw costs if requested
    if store_costs:
        if fmt in ("text", "both"):
            filepath = data_dir / (basename + "_costs.txt")
            write_costs_text(filepath, costs_list, args_dict)
            print(f"Wrote costs text file: {filepath}")

        if fmt in ("hdf5", "both"):
            filepath = data_dir / (basename + "_costs.h5")
            write_costs_hdf5(filepath, costs_list, args_dict, seeds)
            print(f"Wrote costs HDF5 file: {filepath}")

    # Write averaged N(c',d,c) distribution if computed
    if compute_homodist:
        filepath = data_dir / (basename + "_homodist.h5")
        write_homodist_hdf5(filepath, homodist_accumulator.mean, args_dict, seeds, len(seeds))
        print(f"Wrote homogeneous distribution HDF5 file: {filepath}")
        print(f"  Shape: {homodist_accumulator.mean.shape} (num_costs, num_distances, num_costs)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Collects the number of occurences of each cost N(c) for Watts-Strogatz small-world graphs.",
    )

    parser.add_argument("-n", "--numNodes", type=int, required=True, help="Number of nodes/vertices in the graph.")
    parser.add_argument("-k", "--nearestNeighbors", type=int, required=True, help="Each node is joined with its k nearest neighbors in a ring topology.")
    parser.add_argument("-p", "--rewiringProbability", type=float, required=True, help="Probability of rewiring each edge.")
    parser.add_argument("-s", "--seedStart", type=int, required=True, help="Start of the range of seeds.")
    parser.add_argument("-g", "--numGraphs", type=int, required=True, help="Number of graphs to use (seeds will be contiguous range).")
    parser.add_argument("-b", "--backend", default="auto", choices=["auto", "python", "c", "gpu", "gpumpi"], type=str, help="Backend to use for computing maxcut costs. Use 'gpu' for GPU-accelerated N(c',d,c) computation.")
    parser.add_argument("-f", "--format", default="hdf5", choices=["text", "hdf5", "both"], type=str, help="Output format: text (human-readable), hdf5 (compressed, fast), or both.")
    parser.add_argument("-c", "--storeCosts", action="store_true", help="Additionally store raw costs (cost of each bitstring) in separate file(s).")
    parser.add_argument("-d", "--computeHomodist", action="store_true", help="Compute averaged N(c',d,c) homogeneous distribution using Julia. Uses GPU if backend is 'gpu'.")

    try:
        args = parser.parse_args()
        main(args)

    except SystemExit as e:
        print("\nArgument parser failed. Did you provide the correct args?")
        parser.print_help()
