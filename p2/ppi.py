"""
ppi.py
This file uses the node2vec algorithm to embed gene nodes and find similar genes. 
Author: Christopher Sun

Usage: python ppi.py diseaseGeneFile interactionNetworkFile
"""

import sys
import networkx as nx
from node2vec import Node2Vec as n2v
import numpy as np

class PPI(object):
    """
    Object to embed gene nodes and calculate their similarity with a list of diseased genes
    """

    # Function to load the disease gene file and the interaction network file
    # Inputs:
    #       diseaseGeneFile: path to disease gene file
    #       interactionNetworkFile: path to interaction network file
    # Returns:
    #       None
    def load_data(self, diseaseGeneFile, interactionNetworkFile):
        self.disease_genes = []
        with open(diseaseGeneFile, "r") as f:
            for line in f:
                self.disease_genes.append(line.strip())
        self.interactions = []
        with open(interactionNetworkFile, "r") as f:
            for line in f:
                contents = line.split()
                self.interactions.append((contents[0], contents[1]))
        G_interactions = nx.Graph()
        G_interactions.add_edges_from(self.interactions)
        self.interactions = G_interactions
    
    # Function to calculate the embedding of each node using node2vec
    # Inputs:
    #       None
    # Returns:
    #       a tuple containing the list of nodes and their corresponding embeddings
    def calculate_embedding(self):
        g_emb = n2v(self.interactions, dimensions=64, walk_length=30, num_walks=100, workers=1, seed=42)
        model = g_emb.fit(window=3, min_count=1, batch_words=4, workers=1, seed=42)
        embeddings = [model.wv.get_vector(str(node)) for node in self.interactions.nodes]
        return list(self.interactions.nodes), embeddings

    # Function to get a set of genes that are close to disease genes based on a threshold
    # Inputs:
    #       gene_nodes: list of nodes
    #       gene_embeddings: list of embeddings corresponding to gene_nodes
    #       threshold: upper threshold to be considered similar
    # Returns:
    #       a set of similar genes
    def get_close_genes(self, gene_nodes, gene_embeddings, threshold):
        similar_genes = set(self.disease_genes)
        disease_gene_embeddings = np.array([gene_embeddings[i] for i in range(len(gene_embeddings)) if gene_nodes[i] in self.disease_genes])
        gene_embeddings_norm = gene_embeddings / np.linalg.norm(gene_embeddings, keepdims=True, axis=1)
        disease_gene_embeddings_norm = disease_gene_embeddings / np.linalg.norm(disease_gene_embeddings, keepdims=True, axis=1)
        distances_all = 1 - np.dot(gene_embeddings_norm, disease_gene_embeddings_norm.T)
        min_distances = np.min(distances_all, axis=1)
        for i in range(len(min_distances)):
            if min_distances[i] <= threshold:
                similar_genes.add(gene_nodes[i])
        return similar_genes

def main():
    diseaseGeneFile = sys.argv[1]
    interactionNetworkFile = sys.argv[2]
    ppi = PPI()
    ppi.load_data(diseaseGeneFile, interactionNetworkFile)
    genes_all_thresholds = []
    for threshold in [0.05, 0.1, 0.2, 0.3]:
        nodes, embeddings = ppi.calculate_embedding()
        genes_all_thresholds.append(ppi.get_close_genes(nodes, embeddings, threshold))
    print([len(i) for i in genes_all_thresholds])
    
if __name__ == "__main__":
    main()