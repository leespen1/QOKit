import grips, numpy as np, networkx as nx, argparse

parser = argparse.ArgumentParser(
    description="Collects the homogeneous distribution"
)
parser.add_argument("-t", "--type", type=str, default="ErdosRenyi", help="Type of graph to generate (ErdosRenyi)")





