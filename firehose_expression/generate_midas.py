#!/usr/bin/env python

import sys
sys.path.insert(0, "../firehose_mutation")

import csv
import numpy as np
from numpy import *
import firehose_mutation.get_mutation_data as m

samples = []
sample_vals = []
sample_genes = {}

(all_mutations, patient_mutations, patients) = m.get_mutation_data("/Users/sagrava/Work/bioinformatics/ctdd/2_01_2013_deadline/firehose_mutation")

print all_mutations

"""
# get the genes in haian's list that are also in the pathway
gene_list = []
with open('../in_p53_and_haian.txt', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    for row in spamreader:
        gene_list.append(row[0])

# get all the samples from the dataset @samples, 
# the gene and the readout values @sample_genes
with open('LUAD.txt', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    header = spamreader.next()
    header = header[1:]
    j = 0
    for sample in header:
        if j % 3 == 0:
            samples.append(sample)
        j = j + 1

    spamreader.next()
    i = 0
    for row in spamreader:
        #print ', '.join(row)
        (gene, _) = row[0].split("|")
        if gene in gene_list:

            remaining = row[1:]
            j = 1
            for val in remaining:
                if j %3 == 0:
                    sample_genes.setdefault(gene,[]).append(val)
                    #sample_vals.append(val)
                j = j + 1


# initialize the matrix that holds the RPKM values
# multiple by two since we have two time points
# (time point 1 will be initial state of 0 and time point 2
#   will be the actual gene expression values)
data = np.zeros((len(samples)*2, len(gene_list)*2))
data_val_index = len(gene_list) 

time_point_2 = len(samples) 

print data.shape, time_point_2, data_val_index
data[time_point_2:,0:data_val_index] = 30

j = len(gene_list)
for gene in gene_list:
    vals = sample_genes[gene]

    # start i at the next time point
    # everything before that is initialize at 0
    i = len(samples)
    for rpkm in vals:
        data[i][j] = rpkm
        i = i + 1
    j = j + 1


np.savetxt("foo.csv", data, fmt="%.4f", delimiter="\t")

midas_data = []
with open('foo.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    for row in spamreader:
        midas_data.append(row)

header = []
row = []

#header.append("TR:Mock:CellType")
for gene in gene_list:
    row.append("DA:%s" %(gene))

for gene in gene_list:
    row.append("DV:%s" %(gene))

header.append(row)
header.extend(midas_data)

import itertools

myfile = csv.writer(open("midas.csv", "wb"))
for row in header:
    myfile.writerow(row)
"""
