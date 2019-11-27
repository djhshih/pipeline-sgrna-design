#!/usr/bin/env python3

import argparse

pr = argparse.ArgumentParser("Convert oligo tsv file to fasta")
pr.add_argument("input", help="input tsv file [name, sequence]")
pr.add_argument("output", help="output fasta file")

argv = pr.parse_args()

infname = argv.input
outfname = argv.output

outf = open(outfname, 'w')

with open(infname, 'r') as inf:
    header = inf.readline()
    for line in inf:
        tokens = line.split('\t')
        outf.write(">{}\n{}\n".format(tokens[0], tokens[1]))

outf.close()
