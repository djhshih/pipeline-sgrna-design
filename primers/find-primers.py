#!/usr/bin/env python3

class Region:
    '''Genomic region'''

    def __init__(self, chrom=None, start=None, end=None):
        self.chrom = chrom
        self.start = start
        self.end = end

    def from_str(self, region):
        [self.chrom, rest] = region.split(":")
        [start, end] = rest.split("-")
        self.start = int(start)
        self.end = int(end)

    def pad(self, upstream, downstream):
        if self.start > upstream:
            self.start -= upstream
        else:
            self.start = 1
        # NB  end may exceed chromosome length
        self.end += downstream

    def __len__(self):
        return self.end - self.start + 1

    def __str__(self):
        return "{}:{}-{}".format(self.chrom, str(self.start), str(self.end))


class PrimerConfig:
    '''Configuration for Primer3'''

    def __init__(self, id, template, sequence_target):
        self.id = id
        self.template = template
        self.sequence_target = sequence_target
    
    def __str__(self):
        return "\n".join([
            "SEQUENCE_ID={}".format(self.id),
            "SEQUENCE_TEMPLATE={}".format(self.template),
            "PRIMER_PICK_LEFT_PRIMER=1",
            "PRIMER_PICK_RIGHT_PRIMER=1",
            "PRIMER_PRODUCT_SIZE_RANGE=150-900",
            "PRIMER_EXPLAIN_FLAG=0",
            "PRIMER_NUM_RETURN=10",
            "SEQUENCE_TARGET={},{}".format(
                self.sequence_target[0], self.sequence_target[1]),
            "="
            ])

class Primer:
    def __init__(self, sequence, tm):
        self.sequence = sequence
        self.tm = tm

    def __str__(self):
        return self.sequence

class PrimerPair:
    def __init__(self, left, right, product_size):
        self.left = left
        self.right = right
        self.product_size = product_size

    def serialize(self, sep='\t'):
        return '\n'.join([
            sep.join([left.sequence, left.tm]),
            sep.join([right.sequence, right.tm]) ])
 

class Primers:
    '''Container for Primer3 output'''

    def __init__(self, s):
        lines = s.split('\n')
        self.id = self.get_field(lines, "SEQUENCE_ID")
        n = int(self.get_field(lines, "PRIMER_PAIR_NUM_RETURNED"))
        self.n = n
        pairs = []
        for i in range(n):
            pairs.append( self.get_primer_pair(lines, i) )
        self.pairs = pairs

    def primers_iter(self):
        i = 0
        for pair in self.pairs:
            yield (self.id + "-F" + str(i), pair.left.sequence, pair.left.tm)
            yield (self.id + "-R" + str(i), pair.right.sequence, pair.right.tm)
            i += 1

    def pairs_iter(self):
        i = 0
        for pair in self.pairs:
            yield (self.id + "-" + str(i), pair.product_size)
            i += 1

    def get_primer_pair(self, lines, i):
        for line in self.filter(lines, "_{}_".format(i)):
            if not line: continue
            [k, v] = line.split('=')
            if k.find("PRIMER_LEFT_") >= 0:
                if k.find("_SEQUENCE") >= 0:
                    left_sequence = v
                elif k.find("_TM") >=0 :
                    left_tm = v
            elif k.find("PRIMER_RIGHT_") >= 0:
                if k.find("_SEQUENCE") >= 0:
                    right_sequence = v
                elif k.find("_TM") >= 0:
                    right_tm = v
            elif k.find("PRIMER_PAIR_") >= 0:
                if k.find("_PRODUCT_SIZE") >= 0:
                    product_size = v
        return PrimerPair(
            Primer(left_sequence, left_tm),
            Primer(right_sequence, right_tm),
            product_size)

    def filter(self, lines, substr):
        for line in lines:
            if line.find(substr) >= 0:
                yield line

    def get_field(self, lines, key):
        for line in lines:
            [k, v] = line.split('=')
            if k == key:
                return v
        return None

    def __len__(self):
        return self.n


def rm_whitespace(x):
    import string
    return x.translate(str.maketrans('', '', string.whitespace))


import argparse
import subprocess as sp

pr = argparse.ArgumentParser("Find primers to flank target sequence")
pr.add_argument("input", help="input tsv file (name and region coordinate)")
pr.add_argument("ref", help="reference sequence fasta file")
pr.add_argument("primers", help="output primer oligonucleotides tsv file")
pr.add_argument("pairs", help="output primer pairs tsv file")
pr.add_argument("--primer3-settings", help="primer 3 settings file", default="primer3_settings.txt")

argv = pr.parse_args()

infname = argv.input
reffname = argv.ref

primer_padding = 500
target_padding = 50

primersf = open(argv.primers, 'w')
primersf.write("name\tsequence\ttm\n")

pairsf = open(argv.pairs, 'w')
pairsf.write("name\tproduct_size\ttarget_name\n")

with open(infname) as inf:
    header = inf.readline()
    for line in inf:
        [target_name, region_str, strand] = line.rstrip().split('\t')
        pair_name = target_name.replace("-sg", "-e")
        region = Region()
        region.from_str(region_str)
        target_len = len(region)
        region.pad(primer_padding, primer_padding)
        sequence_target = (primer_padding - target_padding + 1, target_len + 2*target_padding)
        out = sp.run(["samtools", "faidx", reffname, str(region)], stdout=sp.PIPE).stdout.decode()
        # omit first line (header) and remove whitespace
        seq = rm_whitespace( out[(out.find('\n')+1):] )
        primer_config = PrimerConfig(pair_name, seq, sequence_target)

        primer_out = sp.run(
            ["primer3_core", "--p3_settings_file={}".format(argv.primer3_settings)],
            input=str(primer_config).encode(), stdout=sp.PIPE
            ).stdout.decode()
        
        primers = Primers(primer_out)
        for p in primers.primers_iter():
            primersf.write('\t'.join(p) + '\n')
        for p in primers.pairs_iter():
            pairsf.write('\t'.join(p) + '\t' + target_name + '\n')

primersf.close()
pairsf.close()

