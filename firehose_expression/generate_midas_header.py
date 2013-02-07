#!/usr/bin/env python

import csv
import numpy as np
from numpy import *

samples = []
sample_vals = []
sample_genes = {}

gene_list = []
with open('../p53.txt', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    for row in spamreader:
        gene_list.append(row[0])

with open('LUAD.txt', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    header = spamreader.next()
    header = header[1:-1]
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

            remaining = row[1:-1]
            j = 1
            for val in remaining:
                if j %3 == 0:
                    sample_genes.setdefault(gene,[]).append(val)
                    #sample_vals.append(val)
                j = j + 1


data = np.zeros((len(samples), len(gene_list)))
data_aq = np.zeros((len(samples), len(gene_list)), dtype=np.uint8)

i = 0

print ("TRT:NONE"),
for gene in gene_list:
    print ("\tDA:%s" %(gene)),

for gene in gene_list:
    #print ("\tDA:%s\tDV:%s" %(gene,gene)),
    print ("DV:%s\t" %(gene)),

print ("\n")
j = 0
for gene in gene_list:
    vals = sample_genes[gene]
    i = 0
    for rpkm in vals:
        data[i][j] = rpkm
        i = i + 1
    j = j + 1

#np.savetxt("foo.csv", data, delimiter="\t")
np.savetxt("zeros.csv", data_aq, fmt='%d',delimiter="\t")
