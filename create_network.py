#!/usr/bin/env python

import csv
import sys

gene_list = []
with open(sys.argv[1]) as f:
    gene_list = f.read().splitlines()

f = open('homo-sapiens-9606.sif', 'rt')
network = []
try:
    reader = csv.reader(f)
    reader.next()
    for row in reader:
        (gene1, interaction, gene2) = row[0].split("\t")
        if interaction == "INTERACTS_WITH" and (gene1 in gene_list and gene2 in gene_list):
            network.append((gene1,gene2))
finally:
    f.close()



for (gene1, gene2) in network:
    print "%s\t1\t%s" % (gene1, gene2)
