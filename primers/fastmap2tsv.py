#!/usr/bin/env python3

import argparse

pr = argparse.ArgumentParser("Convert from fastmap format to tsv")
pr.add_argument("input", help="input fastmap file")
pr.add_argument("output", help="output tsv file")

argv = pr.parse_args()

infname = argv.input
outfname = argv.output

outf = open(outfname, 'w')
outf.write("name\tregion\tstrand\tqlen\tmstart\tmend\tnmatches\tmatches\n")

with open(infname, 'r') as inf:
    while True:
        sq = inf.readline().rstrip()
        if not sq: break
        sq_tokens = sq.split('\t')
        name = sq_tokens[1]
        qlen = sq_tokens[2]

        em = inf.readline().rstrip()
        if not em: break
        em_tokens = em.split('\t')
        mstart = int(em_tokens[1])
        mend = int(em_tokens[2])
        mlen = mend - mstart 
        nmatches = int(em_tokens[3])
        matches = em_tokens[4:]

        # parse the first region
        region = em_tokens[4]
        if not region == '*':
            rtokens = region.split(':')
            chrom = rtokens[0]
            rpos = rtokens[1]
            strand = rpos[0]
            start = int(rpos[1:])
            end = start + mlen - 1
            region = "{}:{}-{}".format(chrom, start, end)
        else:
            region = ""
            strand = ""
            # consume extra blank line that follows '*'
            inf.readline()

        outf.write( "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            name, region, strand, qlen, mstart+1, mend, nmatches,
            ','.join(matches) ) )

        sep = inf.readline().rstrip()
        if not sep: break
        if not sep == "//":
            raise Exception("Missing expected separator '//'")
        
outf.close()

