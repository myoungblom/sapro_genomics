#!/usr/bin/env/python3

import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

####
# Converts core genome alignment into XMFA using a gff file.
####

if len(sys.argv) != 3:
    print("Usage: cfmlCoreAln.py core_alignment_header.gff core_gene_alignment.aln")
    sys.exit(0)

gff = sys.argv[1]
aln = sys.argv[2]

coords = {}

with open(gff,"r") as f:
    for line in f:
        if not line.startswith("#"):
            info = line.strip().split("\t")
            gene = info[8].split(";")[0].split("=")[1]
            start = int(info[3])
            end = int(info[4])
            coords[gene] = [start,end]

seqs = defaultdict(list)

for core in SeqIO.parse(aln,"fasta"):
    for gene,coord in coords.items():
        gene_seq = SeqRecord(Seq(str(core.seq[coord[0]-1:coord[1]])),id=core.id,\
                description="")
        seqs[gene].append(gene_seq)

xmfa = open("clonalframeml.xmfa","w")

for gene,seqs in seqs.items():
    aln = MultipleSeqAlignment(seqs)
    SeqIO.write(aln,xmfa,"fasta")
    xmfa.write("="+gene+"\n")

