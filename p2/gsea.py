"""
gsea.py
This file implements the Gene Set Enrichment Analysis algorithm. 
Author: Christopher Sun

Usage: python gsea.py expfile sampfile keggfile
"""

import sys
import pandas as pd
import numpy as np

class GSEA(object):
    """
    Object to run gene set enrichment analysis on a data set of expression values (# of genes by # of samples) and various gene sets
    """

    # Function to load the expression file, sample files, and gene set file
    # Inputs:
    #       expfile: path to expression file
    #       sampfile: path to sample file
    #       genesets: path to gene set file
    # Returns:
    #       None
    def load_data(self, expfile, sampfile, genesets):
        self.exp_data = pd.read_csv(expfile, sep="\t", index_col=0)
        self.samp_data = {}
        with open(sampfile, "r") as f:
            for line in f:
                contents = line.split()
                self.samp_data[contents[0]] = int(contents[1])
        self.geneset_data = {}
        with open(genesets, "r") as f:
            for line in f:
                contents = line.split()
                self.geneset_data[contents[0]] = contents[2:]

    # Function to rank the genes by level of differential expression (highest to lowest)
    # Inputs:
    #       None
    # Returns:
    #       list of genes ranked in descending order according to log fold change
    def get_gene_rank_order(self):
        healthy_cols = []
        sick_cols = []
        for i in range(len(self.exp_data.columns)):
            if self.samp_data[self.exp_data.columns[i]] == 1:
                sick_cols.append(self.exp_data.columns[i])
            else:
                healthy_cols.append(self.exp_data.columns[i])
        logFC = self.exp_data[sick_cols].mean(axis=1) - self.exp_data[healthy_cols].mean(axis=1)
        return list(logFC.sort_values(ascending=False).index), logFC
    
    # Function to calculate the enrichment score (ES) of a certain geneset
    # Inputs:
    #       geneset: the name of the gene set to calcualate the ES for
    # Returns:
    #       ES rounded to two decimal places
    def get_enrichment_score(self, geneset):
        running_sum = 0
        ES = 0
        gene_rank, _ = self.get_gene_rank_order()
        geneset = [i for i in self.geneset_data[geneset] if i in gene_rank]
        N = len(gene_rank)
        G = len(geneset)
        for i in range(len(gene_rank)):
            if gene_rank[i] in geneset:
                running_sum += ((N - G) / G) ** 0.5
            else:
                running_sum -= (G / (N - G)) ** 0.5
            if running_sum > ES:
                ES = running_sum
        return round(ES, 2)
    
    # Function that runs 100 trials of ES calculations after random permutations of labels
    # and calculates adjusted p-values of each gene set
    # Inputs:
    #       p: threshold p-value before Bonferroni correction is applied
    # Returns:
    #       list of gene sets that are significant (enriched)
    def get_sig_sets(self, p):
        num_trials = 100
        actual_ES = {key: self.get_enrichment_score(key) for key, value in self.geneset_data.items()}
        counts = {key: 0 for key, value in self.geneset_data.items()}
        samp_data_copy = self.samp_data
        for i in range(num_trials):
            labels = np.random.randint(2, size=len(self.samp_data))
            self.samp_data = {key: labels[i] for i, key in enumerate(self.samp_data)}
            for key, value in self.geneset_data.items():
                counts[key] += (self.get_enrichment_score(key) >= actual_ES[key])
        for key, value in counts.items():
            counts[key] /= num_trials
        self.samp_data = samp_data_copy
        return [key for key, value in counts.items() if value <= p / len(self.geneset_data)]

def main():
    expfile = sys.argv[1]
    sampfile = sys.argv[2]
    genesets = sys.argv[3]
    gsea = GSEA()
    gsea.load_data(expfile, sampfile, genesets)
    logfc_idxs, logfc = gsea.get_gene_rank_order()
    print(gsea.get_enrichment_score(['ADRB2',
 'AKT1',
 'AR',
 'AVP',
 'C3',
 'CASP3',
 'CTNNB1',
 'IL13',
 'IL1B',
 'IL6',
 'KRAS',
 'MAPK1',
 'MAPK14',
 'MAPK3',
 'MMP9',
 'NFKB1',
 'POMC',
 'PTGER1',
 'RELA',
 'TGFB1',
 'TIMP1',
 'TP53']))
    # print(logfc["BMP4"])

if __name__ == "__main__":
    main()