"""
networkgen.py
This file generates the edges of a graph of ligands that have statistically significant Tanimoto scores.
Author: Christopher Sun

Usage: python networkgen.py drugs_path targets_path protein_nodes_path
"""

import sys
from chemoUtils import *

# Function that calls bootstrap() to find statistically significant pairs of proteins
def main():
    drugs_path = sys.argv[1]
    targets_path = sys.argv[2]
    protein_nodes_path = sys.argv[3]
    drugs = load_drugs(drugs_path)
    ligands = load_ligands(targets_path)
    protein_nodes = load_protein_nodes(protein_nodes_path)
    tanimoto_scores = load_tanimoto()
    pairs = set()
    with open("network_edgelist.txt", "w") as f:
        for protein_a in protein_nodes:
            for protein_b in protein_nodes:
                if protein_a != protein_b:
                    protein_a_sorted, protein_b_sorted = sorted((protein_a, protein_b))
                    if not (protein_a_sorted, protein_b_sorted) in pairs:
                        pairs.add((protein_a_sorted, protein_b_sorted))
                        p_value = bootstrap(protein_a_sorted, protein_b_sorted, drugs, ligands, tanimoto_scores, 500)
                        if p_value <= 0.05:
                            f.write(f"{protein_a_sorted} {protein_b_sorted}\n")

if __name__ == "__main__":
    main()