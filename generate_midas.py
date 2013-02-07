#!/usr/bin/env python

import sys


import csv
import numpy as np
from numpy import *
from firehose_mutation import get_mutation_data

def get_tumor_suppressor_genes():
    gene_list = []
    with open('tumor_suppressor_list.txt', 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')
        for row in spamreader:
            gene_list.append(row[1])

    return gene_list

def get_oncogenes():
    gene_list = []
    with open('oncogene_list.txt', 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')
        for row in spamreader:
            gene_list.append(row[1])

    return gene_list

def get_treatment_label(gene,i):
    if i == -1:
        return gene+"i"
    elif i == 1:
        return gene
    return "0"

oncogenes = get_oncogenes()
tspgenes  = get_tumor_suppressor_genes()

def calc_treatment(mutations):
    """
        type's of mutations:
            5'UTR
            Frame_Shift_Del
            Frame_Shift_Ins
            Intron
            Missense_Mutation
            Nonsense_Mutation
            Silent
            Splice_Site
    """

    treatments = {}
    treatment = 0
    for (gene, clazz, type) in mutations:
        if clazz== "5'UTR":
            treatment = 0
        elif clazz== "Frame_Shift_Del":
           treatment = -1
        elif clazz== "Frame_Shift_Ins":
           treatment = 1
        elif clazz== "Intron":
           treatment = 0
        elif clazz== "Missense_Mutation":
            if gene in oncogenes:
                treatment = 1 
            elif gene in tspgenes:
                treatment = -1 
        elif clazz== "Nonsense_Mutation":
           treatment = -1 
        elif clazz== "Silent":
            treatment = 0
        elif clazz== "Splice_Site":
            if gene in oncogenes:
                treatment = 1 
            elif gene in tspgenes:
                treatment = -1 
        else:
            treatment = 0

        treatments[gene] = treatment

    return treatments


samples = []
sample_vals = []
sample_genes = {}

(mutation_histogram, patient_mutations, patients) = \
    get_mutation_data.get_all_mutations \
        ("firehose_mutation")


# get the genes in haian's list that are also in the pathway
gene_list = []
#with open('in_p53_and_haian.txt', 'rb') as csvfile:
with open('common_genes.txt', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    for row in spamreader:
        gene_list.append(row[0])

"""
# get all the samples from the dataset @samples, 
# the gene and the readout values @sample_genes
with open('firehose_expression/LUAD.txt', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    header = spamreader.next()
    header = header[1:]
    j = 0
    identifier = None
    for sample in header:
        if j % 3 == 0:
            identifier = sample.rsplit("-",3)[0]
            samples.append(identifier)
        j = j + 1

    print len(samples)
    spamreader.next()
    for row in spamreader:
        (gene, _) = row[0].split("|")
        if gene in gene_list:
            remaining = row[1:]
            i = 0
            j = 1
            for val in remaining:
                if j %3 == 0:
                    identifier = samples[i]
                    sample_genes.setdefault(gene,[]).append((val, identifier))
                    i = i + 1
                j = j + 1

# get patient treatment types
patient_treatments = {}
for id in samples:

    if id in patients:
        #print id, "\n", 
        treatments = calc_treatment(patient_mutations[id])

    else:
        treatments = {}


    patient_treatments[id] = treatments

    
all_classes = []
treatment_index = {}
index = 0
for id in patient_treatments:
    for gene in patient_treatments[id]:
        if patient_treatments[id][gene] == 1 and gene not in all_classes:
            all_classes.append(gene)
            treatment_index[gene] = index
            index = index + 1
        elif patient_treatments[id][gene] == -1 and gene+"i" not in all_classes:
            all_classes.append(gene+"i")
            treatment_index[gene+"i"] = index
            index = index + 1

print all_classes
all_classes = ["TR:" + val for val in all_classes]

# initialize the matrix that holds the RPKM values
# multiple by two since we have two time points
# (time point 1 will be initial state of 0 and time point 2
#   will be the actual gene expression values)

num_treatments = len(all_classes)
data = np.zeros((len(samples)*2, (len(gene_list)*2)+num_treatments+1) )

time_point_2 = len(samples) 

#print data.shape, time_point_2, data_val_index
data[time_point_2:,num_treatments+1:] = 30
data[:,0] = 1

j = len(gene_list)+num_treatments+1
for gene in gene_list:
    vals = sample_genes[gene]

    # start i at the next time point
    # everything before that is initialize at 0
    i = len(samples)
    for (rpkm,id) in vals:
        data[i][j] = rpkm
        #print id,gene,patient_treatments[id]
        if gene in patient_treatments[id]:
            label = get_treatment_label(gene,patient_treatments[id][gene])
            if label in treatment_index:
                data[i][treatment_index[label]] = 1

        i = i + 1
    j = j + 1

(rows,cols)=data.shape
i = len(samples)
while i < rows:
    row_sums = data[i][len(gene_list)+num_treatments+1:]
    data[i][len(gene_list) + num_treatments+1:] = row_sums/(row_sums.max()+1)
    i = i + 1

np.savetxt("foo.csv", data, fmt="%.4f", delimiter="\t")

midas_data = []
with open('foo.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    for row in spamreader:
        midas_data.append(row)

header = []
row = []

row.append("TR:Mock:CellType")
row.extend(all_classes)
for gene in gene_list:
    row.append("DA:%s" %(gene))

for gene in gene_list:
    row.append("DV:%s" %(gene))

header.append(row)
header.extend(midas_data)

import itertools

myfile = csv.writer(open(sys.argv[1], "wb"))
for row in header:
    myfile.writerow(row)
"""
