#!/bin/bash

min_match_len=23
ref=~/data/gencode/release-32/GRCh38.primary_assembly.genome.fa

bwa fastmap -l $min_match_len $ref ../crispr-targets_dna-repair-genes_sgrnas.fa >
	sgrna-targets_fastmap.txt

./fastmap2tsv.py sgrna-targets_fastmap.txt sgrna-targets_fastmap.tsv

./find-primers.py sgrna-targets_fastmap.tsv $ref sgrna-primers_oligos.tsv sgrna-primers_pairs.tsv

./tsv2fasta.py sgrna-primers_oligos.tsv sgrna-primers_oligos.fa

bwa fastmap -l 18 $ref sgrna-primers_oligos.fa \
	> sgrna-primers_fastmap.txt

./fastmap2tsv.py sgrna-primers_fastmap.txt sgrna-primers_fastmap.tsv

./filter-primers.py sgrna-primers_fastmap.tsv sgrna-primers_oligos.tsv sgrna-primers-filtered_oligos.csv

head -n 1 sgrna-primers-filtered_oligos.csv > sgrna-primers-picked_oligos.csv
grep -E '1.?$' sgrna-primers-filtered_oligos.csv >> sgrna-primers-picked_oligos.csv
