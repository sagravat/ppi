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

    for (gene, clazz, type) in mutations:
        if type == "5'UTR":
            treatment = 0
        elif type == 'Frame_Shift_Del':
           treament = -1
        elif type == 'Frame_Shift_Ins':
           treament = 1
        elif type == 'Intron':
           treament = 0
        elif type == 'Missense_Mutation':
            if gene in oncogenes:
                treament = 1 
            elif gene in tspgenes:
                treament = -1 
        elif type == 'Nonsense_Mutation':
           treament = -1 
        elif type == 'Silent':
            treatment = 0
        elif type == 'Splice_Site':
            if gene in oncogenes:
                treament = 1 
            elif gene in tspgenes:
                treament = -1 


        else:
            print "Neither: ", gene, clazz, type


samples = []
sample_vals = []
sample_genes = {}

(mutation_histogram, patient_mutations, patients) = \
    get_mutation_data.get_all_mutations \
        ("/Users/sagrava/Work/bioinformatics/ctdd/2_01_2013_deadline/firehose_mutation")


# get the genes in haian's list that are also in the pathway
gene_list = []
with open('in_p53_and_haian.txt', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    for row in spamreader:
        gene_list.append(row[0])

# get all the samples from the dataset @samples, 
# the gene and the readout values @sample_genes
with open('firehose_expression/LUAD.txt', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    header = spamreader.next()
    header = header[1:]
    j = 0
    for sample in header:
        if j % 3 == 0:
            identifier = sample.rsplit("-",3)[0]
            samples.append(identifier)
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
                j = j + 1


# initialize the matrix that holds the RPKM values
# multiple by two since we have two time points
# (time point 1 will be initial state of 0 and time point 2
#   will be the actual gene expression values)
data = np.zeros((len(samples)*2, len(gene_list)*2))
data_val_index = len(gene_list) 

time_point_2 = len(samples) 

#print data.shape, time_point_2, data_val_index
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

(rows,cols) = data.shape
for i in range(rows/2):
    id = samples[i]

    if id in patients:
        #print id, "\n", 
        calc_treatment(patient_mutations[id])
    #else:
        #print id, [],


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
