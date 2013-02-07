#!/usr/bin/env python

import csv
import glob
import os


# get the genes in haian's list that are also in the pathway
def get_all_mutations(dir):
    gene_list = []
    with open('common_genes.txt', 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')
        for row in spamreader:
            gene_list.append(row[0])
    
    patients = []
    patient_mutations = {}
    patient_mutation_genes = {}
    # histogram
    mutation_histogram = {}
    patient_id = None
    mutated_genes = []
    for files in glob.glob("%s/*.txt" %dir):
        with open(files, 'rb') as csvfile:
            spamreader = csv.reader(csvfile, delimiter='\t')
            spamreader.next()
            for row in spamreader:
                if row[0] in gene_list:
                    gene = row[0]
                    if gene not in mutated_genes:
                        mutated_genes.append(gene)
                    variant_class = row[8]
                    variant_type  = row[9]
                    patient_id = row[15].rsplit("-",3)[0]
                    patient_mutations.setdefault(patient_id, []).append([gene,variant_class,variant_type])
                    patient_mutation_genes.setdefault(patient_id, []).append([gene])
                    if gene in mutation_histogram:
                        mutation_histogram[gene] =  mutation_histogram[gene] + 1
                    else:
                        mutation_histogram[gene] = 1
            patients.append(patient_id)

    for gene in mutated_genes:
        print gene
    print patient_mutation_genes
    return (mutation_histogram,patient_mutations, patients)


