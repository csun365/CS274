"""
pvalue.py
This file calculates the significance of the Tanimoto score given two proteins.
Author: Christopher Sun

Usage: python pvalue.py -n num_iterations -r seed drugs_path targets_path protein_a protein_b
"""

import sys
from chemoUtils import *
from tanimoto import *

# Function that calls bootstrap() to calculate the p-value of the Tanimoto score between two proteins
def main():
    args = sys.argv[1:]
    num_iterations = 500
    seed = 214
    if "-n" in args:
        num_iterations = int(args[args.index("-n") + 1])
    if "-r" in args:
        seed = int(args[args.index("-r") + 1])
    drugs_path = args[-4]
    targets_path = args[-3]
    protein_a = args[-2]
    protein_b = args[-1]
    drugs = load_drugs(drugs_path)
    ligands = load_ligands(targets_path)
    # tanimoto_scores = load_tanimoto(generate_all_tanimoto(drugs, ligands, output_path=None))
    tanimoto_scores = load_tanimoto()
    p_value = bootstrap(protein_a, protein_b, drugs, ligands, tanimoto_scores, num_iterations, seed)
    print(p_value)
    
if __name__ == "__main__":
    main()