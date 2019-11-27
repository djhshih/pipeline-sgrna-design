#!/usr/bin/env python3

import csv
import sys

from Bio.Seq import Seq

min_target_len = 17

infname = "crispr-targets_dna-repair-genes_sgrnas.tsv"
outfname = "crispr-targets_dna-repair-genes_oligos.csv"

outf = open(outfname, 'w', newline='')
writer = csv.DictWriter(outf, fieldnames=["name", "sequence"], delimiter=',')
writer.writeheader()

with open(infname, newline='') as inf:
    reader = csv.DictReader(inf, delimiter='\t')
    for r in reader:
        target = r["target_sequence"]

        # do not create oligos for sgRNA with plasmid already available
        if len(r["name"]) > 0: continue

        # trim 5' of target as necessary to ensure first base is a G,
        # as required by the U6 promoter
        target = target[target.find('G'):]
        if len(target) < min_target_len:
            print("{} is too short after trimming\n".format(r["oligo_name"]), file = sys.stderr)
            continue

        target_seq = Seq(target)
        target_seq_rc = target_seq.reverse_complement()
        oligo1 = Seq("CACC") + target_seq
        oligo2 = Seq("AAAC") + target_seq_rc
        pattern = "GN{}NGG".format(len(target))
        oligo_name = r["oligo_name"]
        
        writer.writerow({"name": oligo_name + "-F", "sequence": str(oligo1)})
        writer.writerow({"name": oligo_name + "-R", "sequence": str(oligo2)})

        print(oligo_name)
        print("5' - " + oligo1 + "     - 3'")
        print("3' -     " + oligo2[::-1] + " - 5'")
        print()
        
outf.close()

