"""
plot_graph.py
This file plots the graph of ligands and colors nodes by indications.
Author: Christopher Sun

Usage: python plot_graph.py network_edgelist_path protein_nodes_path output_path
"""
import sys
import networkx
import matplotlib.pyplot as plt
from chemoUtils import *

# Function that uses networkx to construct and plot network of proteins and indications
def main():
    network_edgelist_path = sys.argv[1]
    protein_nodes_path = sys.argv[2]
    output_path = sys.argv[3]
    protein_dict, indications_dict = load_protein_dict(protein_nodes_path)
    G_network = networkx.read_edgelist(network_edgelist_path)
    G_network = networkx.relabel_nodes(G_network, protein_dict)
    print(G_network)
    color_dict = {"bp": "red", "bp;cholesterol": "green", "bp;cholesterol;diabetes": "blue", "bp;diabetes": "purple"}
    plt.figure(figsize=(8,8))
    networkx.draw_networkx(G_network, node_color=[color_dict[indications_dict[i]] for i in G_network])
    plt.savefig(output_path, dpi=150)
    plt.show()

if __name__ == "__main__":
    main()