"""
tanimoto.py
This file calculates the Tanimoto scores between every pair of proteins and writes them to a file.
Author: Christopher Sun

Usage: python tanimoto.py drugs_path targets_path output_path
"""

import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
from chemoUtils import *

# Function to calculate tanimoto score between two fingerprints
# Inputs:
#       fpt1: fingerprint 1
#       fpt2: fingerprint 2
# Returns:
#       Tanimoto score with 6 decimal places
def calculate_tanimoto(fpt1, fpt2):
    tanimoto = len(set(fpt1).intersection(fpt2)) / len(set(fpt1).union(fpt2))
    return "{:.6f}".format(round(tanimoto, 6))

# Function to determine if two drugs share a target
# Inputs:
#       drug1: first drug
#       drug2: second drug
#       targets: dictionary of targets
# Returns:
#       indicator variable that represents sharing of a target
def shares_target(drug1, drug2, targets):
    if not drug1 in targets or not drug2 in targets:
        return 0
    return int(len(targets[drug1].intersection(targets[drug2])) > 0)

# Function to generate all Tanimoto scores
# Inputs:
#       drugs: drug file
#       targets: targets file
#       output_path: optional output path of generated csv file
# Returns:
#       dataframe of Tanimoto scores
def generate_all_tanimoto(drugs, targets, output_path=None):
    drug1 = []
    drug2 = []
    scores = []
    shares = []
    for i in range(drugs.shape[0]):
        for j in range(i + 1, drugs.shape[0]):
            drug1.append(drugs["db_id"].iloc[i])
            drug2.append(drugs["db_id"].iloc[j])
            scores.append(calculate_tanimoto(drugs["maccs"].iloc[i], drugs["maccs"].iloc[j]))
            shares.append(shares_target(drugs["db_id"].iloc[i], drugs["db_id"].iloc[j], targets))
    output = pd.DataFrame([drug1, drug2, scores, shares]).T
    if output_path:
        output.to_csv(output_path, index=False, header=False)
    return output

# Function to plot histograms of Tanimoto scores
# Inputs:
#       title: title of histogram
#       save_title: title of png of saved histogram
# Returns:
#       None
def plot_hists(title, save_title):
    tanimoto = pd.read_csv("tanimoto_scores.csv")
    if save_title == "all_tanimoto":
        data = tanimoto.iloc[:,2]
    elif save_title == "shared_tanimoto":
        data = tanimoto[tanimoto.iloc[:,-1] == 1].iloc[:,2]
    elif save_title == "notshared_tanimoto":
        data = tanimoto[tanimoto.iloc[:,-1] == 0].iloc[:,2]
    plt.hist(data, bins=1+math.ceil(np.log2(data.shape[0])))
    plt.title(title)
    plt.xlabel("Tanimoto Score")
    plt.ylabel("Number of Drug Pairs")
    plt.tight_layout()
    plt.savefig(save_title + ".png", dpi=250)
    plt.show()

# Function that calls generate_all_tanimoto() and plot_hists()
def main():
    drugs_path = sys.argv[1]
    targets_path = sys.argv[2]
    output_path = sys.argv[3]
    drugs = load_drugs(drugs_path)
    targets = load_targets(targets_path)
    generate_all_tanimoto(drugs, targets, output_path)
    # plot_hists("csun27 All", "all_tanimoto")
    # plot_hists("csun27 Shared", "shared_tanimoto")
    # plot_hists("csun27 Not Shared", "notshared_tanimoto")

if __name__ == "__main__":
    main()