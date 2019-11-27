# CRISPR sgRNA and primer design

sgRNA was designed using the CHOPCHOP web interface and a tab-delimited table
was prepared.

sgRNA are named based on the target exon number of the RefSeq transcript
with the numerically smallest (earliest) accession number.

sgRNA sequences (with PAMs) were mapped to the human genome to ensure unique 
mapping with `bwa fastmap`. Non-unique sgRNAs were re-designed.

PCR primers were designed against each sgRNA site using Primer3 with 
Primer3Web v4.0.0 default settings.
Forward primers always map to the plus stand of the genome and reverse primers
to the minus strand.
PCR primers are named based on the exon targeted by the sgRNA, and they flank
the target sgRNA site.

See `requirements.txt` for prerequisite software installations.

## Remarks

[CHOPCHOP](https://bitbucket.org/valenlab/chopchop) does not currently use
updated settings for Primer3.
[Benchling](https://benchling.com) does not seem to use the updated
settings either. Therefore, we refrain from designing primers using these tools.

