import grips, argparse, networkx as nx, numpy as np
from collections import OrderedDict
from pathlib import Path


def main(args):
    args_dict = OrderedDict(sorted(vars(args).items()))
    args_dict.pop("backend") # Don't keep track of backend
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
    seeds = range(args.seedStart, args.seedStart + args.numGraphs)
    graphs = [nx.erdos_renyi_graph(n, p, seed=seed) for seed in seeds]


    # Directory of this script
    script_dir = Path(__file__).resolve().parent

    # Data directory relative to script
    data_dir = script_dir / f"Data_graphType={graphType}/numNodes={n}"
    data_dir.mkdir(exist_ok=True)

    # Output file path
    filename = grips.args_to_str(args_dict, "_") + ".txt"
    filepath = data_dir / filename

    with open(filepath, "w") as f:
        for graph in graphs:
            # Get cost for each bitstring
            costs = np.rint(grips.get_costs(graph, args.backend)).astype(int)
            N_c = np.bincount(costs)
            np.savetxt(f, N_c[None, :], fmt="%d")

            ## Binary save version
            ## write length
            #np.array([len(N_c)], dtype=np.int64).tofile(f)
            ## write data
            #N_c.tofile(f)


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

    try:
        args = parser.parse_args()
        main(args)

    except SystemExit as e:
        print("\nArgument parser failed. Did you provide the correct args?")
        parser.print_help()

    

