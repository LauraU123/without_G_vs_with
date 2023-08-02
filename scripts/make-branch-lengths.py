import argparse
from Bio import Phylo
import pandas as pd
import json


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="turn newick trees into branch length jsons",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type = str, help='divtree for branch length.', required=True) 
    parser.add_argument('--output-node-data', help="name of the file to write node data to", required=True)
    args = parser.parse_args()

    divtree = Phylo.read(args.tree, 'newick')

    node_data = {}
    for node in divtree.find_clades():
        node_name = node.name or node.confidence
        node_data[node_name] = {}
        if node_name in node_data:
            print(node.branch_length)
            node_data[node_name]["mutation_length"] = node.branch_length


    with open(args.output_node_data, 'w') as fh:
        json.dump({"nodes": node_data}, fh)