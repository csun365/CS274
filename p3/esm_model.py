"""
esm_model.py
This file generates protein embeddings using ESMFold.
Author: Christopher Sun

Usage: python esm_model.py fasta_file output_dir prot_family_file
"""

import sys
from Bio import SeqIO
import torch
import esm 
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import pickle
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

class ESM_model():
    """
    Object to load protein string sequences, generate embeddings, and pool/visualize these embeddings
    """

    # Function to load the fasta file and return a list of (id, sequence) tuples
    # Inputs:
    #       fasta_file: path to fasta file
    # Returns:
    #       list of (id, sequence) tuples
    def load_fasta(self, fasta_file):
        return [(record.id.split("|")[1], record.seq) for record in SeqIO.parse(fasta_file, "fasta")] 

    # Function to load ESM and run inference on protein sequences
    # Inputs:
    #       list_tuples_protId_seq: list of (id, sequence) tuples
    # Returns:
    #       dictionary of id: embedding pairs, where the embedding is 2D: [seq_length, embedding_dim]
    def get_vectors(self, list_tuples_protId_seq):
        device = torch.device("mps")
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        model = model.to(device)
        batch_converter = alphabet.get_batch_converter()
        model.eval()
        result = {}
        for tup in list_tuples_protId_seq:
            _, _, batch_tokens = batch_converter([tup])
            batch_tokens = batch_tokens.to(device)
            with torch.no_grad():
                results = model(batch_tokens, repr_layers=[33], return_contacts=False)
            token_representations = results["representations"][33]
            result[tup[0]] = token_representations[0:1, 1:-1, :].cpu()
        return result

    # Inputs:
    #       dict_protId_embedding: dictionary of id: embedding pairs
    #       mean_max_param: string specifying pooling mode
    # Returns:
    #       dictionary of id: pooled_embedding pairs, where the pooled_embedding is 1D: [1, embedding_dim]
    def pool_representations(self, dict_protId_embedding, mean_max_param):
        if mean_max_param == "mean":
            return {key: np.mean(value.numpy(), axis=1) for key, value in dict_protId_embedding.items()}
        elif mean_max_param == "max":
            return {key: np.max(value.numpy(), axis=1) for key, value in dict_protId_embedding.items()}
    
    # Inputs:
    #       family_data_file: path to protein family file
    #       pooled_vectors_dict: dictionary of id: pooled_embedding pairs
    # Returns:
    #       dataframe of silhouette scores for each family 
    def calculate_silhouette_score(self, family_data_file, pooled_vectors_dict):
        family_lookup = {}
        with open(family_data_file, "r") as f:
            for i in f.readlines()[1:]:
                family_lookup[i[:i.find(",")]] = i[i.find(",") + 1:].strip()
        silhouette_scores = {}
        X = np.array(list(pooled_vectors_dict.values())).reshape(-1,1280)
        df_family = pd.read_csv(family_data_file)
        for family in df_family["family"].unique():
            labels = [family_lookup[id] == family for id, _ in pooled_vectors_dict.items()]
            silhouette_scores[family] = [silhouette_score(X, labels)]
        return pd.DataFrame.from_dict(silhouette_scores, orient="index").reset_index().rename(columns={"index": "Family", 0: "Score"})

# Inputs:
#       data: multidimensional data to pass into TSNE
#       title: title of resulting figure
# Returns:
#       None
def plot_tsne(self, data, title):
    tsne = TSNE(n_components=2, random_state=42)
    transformed = tsne.fit_transform(data)
    plt.scatter(transformed[:,0], transformed[:,1])
    plt.title(title)
    plt.savefig("outputs/" + title + ".png", dpi=250)
    plt.show()

def main():
    fasta_file = sys.argv[1]
    output_dir = sys.argv[2]
    prot_family_file = sys.argv[3]
    esm_model = ESM_model()
    strings = esm_model.load_fasta(fasta_file)
    print(strings)
    embeddings = esm_model.get_vectors(strings)
    # embeddings = pickle.load(open("outputs/embeddings.pkl", "rb"))
    mean_embeddings = esm_model.pool_representations(embeddings, "mean")
    max_embeddings = esm_model.pool_representations(embeddings, "max")
    pickle.dump(embeddings, open("outputs/embeddings.pkl", "wb"))
    pickle.dump(mean_embeddings, open("outputs/meanPooled_embeddings.pkl", "wb"))
    pickle.dump(max_embeddings, open("outputs/maxPooled_embeddings.pkl", "wb"))
    # mean_embeddings = pickle.load(open("outputs/meanPooled_embeddings.pkl", "rb"))
    # max_embeddings = pickle.load(open("outputs/maxPooled_embeddings.pkl", "rb"))
    # plot_tsne(np.array(list(mean_embeddings.values())).reshape(-1,1280), "maxPooled_viz")
    # plot_tsne(np.array(list(max_embeddings.values())).reshape(-1,1280), "meanPooled_viz")
    # mean_silhouette_scores = esm_model.calculate_silhouette_score(prot_family_file, mean_embeddings)
    # max_silhouette_scores = esm_model.calculate_silhouette_score(prot_family_file, max_embeddings)
    # mean_silhouette_scores.to_csv("outputs/meanPooled_silhouette.csv")
    # max_silhouette_scores.to_csv("outputs/maxPooled_silhouette.csv")
    
if __name__ == "__main__":
    main()