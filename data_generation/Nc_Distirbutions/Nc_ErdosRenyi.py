import grips, argparse, networkx as nx, numpy as np
from collections import OrderedDict
from pathlib import Path


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
    """
    import h5py

    # Pad arrays to same length for efficient storage
    max_len = max(len(N_c) for N_c in Nc_list)
    Nc_padded = np.zeros((len(Nc_list), max_len), dtype=np.int64)
    for i, N_c in enumerate(Nc_list):
        Nc_padded[i, :len(N_c)] = N_c

    with h5py.File(filepath, "w") as f:
        # Store data with gzip compression (good balance of speed/size)
        f.create_dataset("Nc", data=Nc_padded, compression="gzip", compression_opts=4)
        f.create_dataset("seeds", data=np.array(seeds))

        # Store parameters as attributes
        for key, val in args_dict.items():
            f.attrs[key] = val


def main(args):
    args_dict = OrderedDict(sorted(vars(args).items()))
    args_dict.pop("backend")  # Don't keep track of backend
    args_dict.pop("format")   # Don't keep track of format
    graphType = "ErdosRenyi"
    args_dict["graphType"] = graphType

    # Rearrange arg order for filename, make easier to parse
    args_dict.move_to_end("edgeProbability", last=False)
    args_dict.move_to_end("numNodes", last=False)
    args_dict.move_to_end("graphType", last=False)
    args_dict.move_to_end("seedStart", last=True)
    args_dict.move_to_end("numGraphs", last=True)

    print("Running script with the following parameters:")
    print("\t", grips.args_to_str(args_dict, "\n\t"))

    n = args.numNodes
    p = args.edgeProbability
    seeds = list(range(args.seedStart, args.seedStart + args.numGraphs))
    graphs = [nx.erdos_renyi_graph(n, p, seed=seed) for seed in seeds]

    # Directory of this script
    script_dir = Path(__file__).resolve().parent

    # Data directory relative to script
    data_dir = script_dir / f"Data_graphType={graphType}/numNodes={n}"
    data_dir.mkdir(exist_ok=True, parents=True)

    # Compute all N(c) distributions
    Nc_list = []
    for graph in graphs:
        costs = np.rint(grips.get_costs(graph, args.backend)).astype(int)
        N_c = np.bincount(costs)
        Nc_list.append(N_c)

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Collects the number of occurences of each cost N(c) for Erdos-Renyi graphs.",
    )

    parser.add_argument("-n", "--numNodes", type=int, required=True, help="Number of nodes/vertices in the graph.")
    parser.add_argument("-p", "--edgeProbability", type=float, required=True, help="Probability of an edge between each pair of vertices")
    parser.add_argument("-s", "--seedStart", type=int, required=True, help="Start of the range of seeds.")
    parser.add_argument("-g", "--numGraphs", type=int, required=True, help="Number of graphs to use (seeds will be contiguous range).")
    parser.add_argument("-b", "--backend", default="auto", choices=["auto", "python", "c", "gpu", "gpumpi"], type=str, help="Backend to use for computing maxcut costs.")
    parser.add_argument("-f", "--format", default="hdf5", choices=["text", "hdf5", "both"], type=str, help="Output format: text (human-readable), hdf5 (compressed, fast), or both.")

    try:
        args = parser.parse_args()
        main(args)

    except SystemExit as e:
        print("\nArgument parser failed. Did you provide the correct args?")
        parser.print_help()

