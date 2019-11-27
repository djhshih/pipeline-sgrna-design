#!/usr/bin/env python3

import csv
import sys

infname = "crispr-targets_dna-repair-genes_sgrnas.tsv"
outfname = "crispr-targets_dna-repair-genes_sgrnas.fa"

outf = open(outfname, 'w')

with open(infname, newline='') as inf:
    reader = csv.DictReader(inf, delimiter='\t')
    for r in reader:
        oligo_name = r["oligo_name"]
        target = r["target_sequence"]
        pam = r["pam"]

        # target sequence and pam have already been reverse complemented
        # s.t. all sequences are on the plus strand frame

        outf.write(">" + oligo_name + "\n")
        outf.write(target + pam + "\n")
        
outf.close()

