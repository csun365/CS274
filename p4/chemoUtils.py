"""
chemoUtils.py
This file contains helper methods for loading files and computing Tanimoto scores that are called by other scripts.
Author: Christopher Sun
"""

import pandas as pd
import numpy as np

# Function to load drug file
# Inputs:
#       drugs_path: path to drug file
# Returns:
#       dataframe of drugs and their fingerprints
def load_drugs(drugs_path):
    drugs = pd.read_csv(drugs_path)
    drugs["maccs"] = drugs["maccs"].apply(lambda x: [int(i) for i in x.split()])
    return drugs

# Function to load targets file
# Inputs:
#       targets_path: path to targets file
# Returns:
#       dictionary with drugs as key and uniprot id as values
def load_targets(targets_path):
    targets_dict = {}
    targets_df = pd.read_csv(targets_path)
    for i in range(targets_df.shape[0]):
        drug = targets_df["db_id"].iloc[i]
        if not drug in targets_dict:
            targets_dict[drug] = set(list(targets_df[targets_df["db_id"] == drug]["uniprot_id"].values))
    return targets_dict

# Function to load ligands file
# Inputs:
#       targets_path: path to targets file
# Returns:
#       dictionary with proteins as key and db id as values
def load_ligands(targets_path):
    ligands_dict = {}
    ligands_df = pd.read_csv(targets_path)
    for i in range(ligands_df.shape[0]):
        protein = ligands_df["uniprot_accession"].iloc[i]
        if not protein in ligands_dict:
            ligands_dict[protein] = set(list(ligands_df[ligands_df["uniprot_accession"] == protein]["db_id"].values))
    return ligands_dict

# Function to load Tanimoto scores
# Inputs:
#       tanimoto_df: dataframe of Tanimoto scores
# Returns:
#       dictionary with tuple of ligands as key and Tanimoto score as values
def load_tanimoto():
    tanimoto_dict = {}

    with open("tanimoto_scores.csv", "r") as f:
        for line in f:
            parts = line.strip().split(",")
            tanimoto_dict[(parts[0], parts[1])] = float(parts[2])
    
    # for _, row in tanimoto_df.iterrows():
    #     tanimoto_dict[(row[0], row[1])] = float(row[2])
    return tanimoto_dict

# Function to load protein nodes
# Inputs:
#       protein_nodes_path: path to protein nodes file
# Returns:
#       dataframe of protein nodes
def load_protein_nodes(protein_nodes_path):
    return pd.read_csv(protein_nodes_path)["uniprot_accession"].values

# Function to load protein and indication dictionaries
# Inputs:
#       protein_nodes_path: path to protein nodes file
# Returns:
#       dictionary mapping uniprot accession to uniprot id, and dictionary mapping uniprot id to indications
def load_protein_dict(protein_nodes_path):
    protein_dict = {}
    indications_dict = {}
    with open(protein_nodes_path, "r") as f:
        for line in f.readlines()[1:]:
            parts = line.split(",")
            protein_dict[parts[0]] = parts[1]
            indications_dict[parts[1]] = parts[2].rstrip()
    return protein_dict, indications_dict

# Function to compute Tanimoto similarity summary statistic
# Inputs:
#       ligands_a: first set of ligands
#       ligands_b: second set of ligands
#       tanimoto_scores: dictionary mapping ligand pairs to Tanimoto scores
# Returns:
#       Tanimoto similarity summary statistic
def compute_similarity(ligands_a, ligands_b, tanimoto_scores):
    summary = 0
    for drug1 in ligands_a:
        for drug2 in ligands_b:
            if drug1 != drug2:
                drug1_sorted, drug2_sorted = sorted((drug1, drug2))
                score = tanimoto_scores[(drug1_sorted, drug2_sorted)]
                summary += score * (score > 0.5)
            else:
                summary += 1
    return summary 

# Function to conduct bootstrapping for Tanimoto summary score statistically significance
# Inputs:
#       protein_a: first protein
#       protein_b: second protein
#       drugs: drugs dataframe
#       ligands: ligands dictionary
#       tanimoto_scores: dictionary mapping ligand pairs to Tanimoto scores
#       num_iterations: bootstrapping iterations
#       seed: random seed
# Returns:
#       p-value
def bootstrap(protein_a, protein_b, drugs, ligands, tanimoto_scores, num_iterations=500, seed=214):
    np.random.seed(seed)
    T_summary = compute_similarity(ligands[protein_a], ligands[protein_b], tanimoto_scores)
    n_a = len(ligands[protein_a])
    n_b = len(ligands[protein_b])
    p_value = 0
    for i in range(num_iterations):
        ligands_a = drugs.iloc[np.random.choice(drugs.shape[0], size=n_a, replace=True).tolist(),0].values
        ligands_b = drugs.iloc[np.random.choice(drugs.shape[0], size=n_b, replace=True).tolist(),0].values
        similarity = compute_similarity(ligands_a, ligands_b, tanimoto_scores)
        if similarity >= T_summary:
            p_value += 1
    return p_value / num_iterations