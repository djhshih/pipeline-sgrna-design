#!/usr/bin/env python3

import argparse
import csv

pr = argparse.ArgumentParser("Filter primers based on alignment hits")
pr.add_argument("fastmap", help="input fastmap tsv file")
pr.add_argument("primers", help="input primer oligo tsv file")
pr.add_argument("output", help="output oligo tsv file")

argv = pr.parse_args()

infname = argv.primers
fastmapfname = argv.fastmap
outfname = argv.output

selected_dict = dict()

filtered = set()

def oligo_name_stem(x):
    i = x.rfind("-F")
    if i < 0:
        i = x.rfind("-R")
    return x[:i]

with open(fastmapfname, 'r', newline='') as fastmapf:
    reader = csv.DictReader(fastmapf, delimiter='\t')
    while True:
        # forward oligo
        try:
            row1 = next(reader)
            row2 = next(reader)
        except:
            break

        if not row1 or not row2: break

        assert(row1['name'].replace("-F", "") == row2['name'].replace("-R", ""))

        # primers should be on opposite strands
        if row1['strand'] == row2['strand']: continue

        # primer should be fully matched
        if row1['qlen'] != row1['mend'] or row1['mstart'] != '1' or row2['qlen'] != row2['mend'] or row2['mstart'] != '1':
                continue
        
        # both primers should be unique
        if int(row1['nmatches']) > 1 or int(row2['nmatches']) > 1: continue

        name1 = row1['name']
        name2 = row2['name']

        stem1 = oligo_name_stem(name1)
        stem2 = oligo_name_stem(name2)

        assert(stem1 == stem2)

        if not stem1 in selected_dict:
            selected_dict[stem1] = ( name1, name2 )

        filtered.add(row1['name']) 
        filtered.add(row2['name']) 

selected = set([ x for pairs in selected_dict.values() for x in pairs ])

with open(infname, 'r', newline='') as inf, open(outfname, 'w', newline='') as outf:
    reader = csv.DictReader(inf, delimiter='\t')
    fieldnames = list(reader.fieldnames)
    fieldnames.append('pick')
    writer = csv.DictWriter(outf, delimiter=',', fieldnames=fieldnames)
    writer.writeheader()

    for row in reader: 
        if row['name'] in filtered:
            if row['name'] in selected:
                row['pick'] = 1
            else:
                row['pick'] = 0
            writer.writerow(row)

