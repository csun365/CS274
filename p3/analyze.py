"""
analyze.py
This file analyzes protein binding site predictions made by ESM.
Author: Christopher Sun

Usage: python analyze.py predicted_results_dir binding_site_file
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import sys

class AnalyzePredictions():
    """
    Object to analyze protein binding site predictions and compute metrics
    """

    # Function to load all predictions and concatenate them to a dataframe
    # Inputs:
    #       predicted_results_dir: path to prediction files
    #       cutoff: threshold for classification of binding site
    # Returns:
    #       dataframe with binding site predictions
    def load_data(self, predicted_results_dir, cutoff):
        all_prots = pd.DataFrame()
        for file in os.listdir(predicted_results_dir):
            df = pd.read_csv(predicted_results_dir + file, sep="\t")
            df = df[df["Prob"] >= cutoff]
            df["protid"] = file[:6]
            all_prots = pd.concat([all_prots, df])
        return all_prots
    
    # Function to load true binding sites
    # Inputs:
    #       protein_bindingSite_file: path to true binding site file
    # Returns:
    #       dataframe with ground truth binding sites
    def get_bindingSite_labels(self, protein_bindingSite_file):
        return pd.read_csv(protein_bindingSite_file)
    
    # Function to calculate accuracy, precision, true positive rate, and false positive rate of predictions
    # Inputs:
    #       predicted_results_dir: path to prediction files
    #       cutoff: threshold for classification of binding site
    #       bindingSite_dataframe: dataframe with ground truth binding sites
    # Returns:
    #       accuracy, precision, true positive rate, and false positive rate of predictions
    def calculate_metrics(self, predicted_results_dir, cutoff, bindingSite_dataframe):
        predictions = self.load_data(predicted_results_dir, cutoff)
        predictions = predictions[predictions["protid"].isin(bindingSite_dataframe["protid"].unique())]
        probs = self.load_data(predicted_results_dir, 0)
        tp, fp, tn, fn = 0, 0, 0, 0
        for protid in predictions["protid"].unique():
            tp_prot, fp_prot, tn_prot, fn_prot = 0, 0, 0, 0
            subset_predictions = set(list(predictions[predictions["protid"] == protid]["Index"]))
            subset_binding = set(list(bindingSite_dataframe[bindingSite_dataframe["protid"] == protid]["binding_site"]))
            for i in subset_predictions:
                if i in subset_binding:
                    tp_prot += 1
                else:
                    fp_prot += 1
            for i in subset_binding:
                if not i in subset_predictions:
                    fn_prot += 1
            for i in set(list(probs[probs["protid"] == protid]["Index"])) - subset_predictions:
                if not i in subset_binding:
                    tn_prot += 1
            tp += tp_prot
            fp += fp_prot
            fn += fn_prot
            tn += tn_prot
        return (tp + tn) / (tp + fp + tn + fn), tp / (tp + fp), tp / (tp + fn), fp / (fp + tn)

# Function to plot histogram of binding site frequencies for each amino acid
# Inputs:
#       predictions: dataframe of binding site predictions
# Returns:
#       None
def plot_hist(predictions):
    counts = {}
    for aa in predictions["AA"].unique():
        counts[aa] = predictions[predictions["AA"] == aa].shape[0]
    plt.bar(counts.keys(), counts.values())
    plt.ylabel("Frequency")
    plt.savefig("outputs/AA_binding_site_freq.png", dpi=250)
    plt.show()

def main():
    predicted_results_dir = sys.argv[1]
    binding_site_file = sys.argv[2]
    analysis = AnalyzePredictions()
    tprs, fprs = [], []
    for cutoff in [i / 10.0 for i in range(11)]:
        predictions = analysis.load_data(predicted_results_dir, cutoff)
        # plot_hist(predictions)
        binding_df = analysis.get_bindingSite_labels(binding_site_file)
        acc, precision, tpr, fpr = analysis.calculate_metrics(predicted_results_dir, cutoff, binding_df)
        tprs.append(tpr)
        fprs.append(fpr)
    plt.plot(tprs, fprs)
    plt.xlabel("TPR")
    plt.ylabel("FPR")
    plt.title("ROC")
    # plt.savefig("outputs/ROC.png", dpi=250)
    plt.show()

if __name__ == "__main__":
    main()