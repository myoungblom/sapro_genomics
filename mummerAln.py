#!/usr/bin/env python3

import sys
import subprocess
from Bio import SeqIO

#######
# Whole genome alignment to a reference sequence using Mummer.
#######

# check for correct arguments
if len(sys.argv) < 3:
	print("Usage: mummerAln.py <reference.fasta> <genome1> ... <genomeN>")
	sys.exit(0)

reffile = sys.argv[1]
fastas = sys.argv[2:]

def runMummer(genome,reference):
    iso = genome.split(".")[0]+"_aln"
    # run mummer
    subprocess.call(["/opt/PepPrograms/mummer-4.0.0beta2/nucmer","--prefix",iso,reference,genome])
    with open(iso+".snps","w") as snpfile:
        subprocess.call(["/opt/PepPrograms/mummer-4.0.0beta2/show-snps","-CHT",iso+".delta"],stdout=snpfile)
    snpfile.close()
    # read reference into dictionary of string(s)
    refDict = {}
    for contig in SeqIO.parse(reference,"fasta"):
        refDict[contig.id] = list(str(contig.seq))
    # read mummer output, change reference sequence accordingly
    with open(iso+".snps","r") as f:
        for line in f:
            info = line.strip().split("\t")
            refcontig = info[8]
            if info[1] != ".":
                position = int(info[0])
                if info[2] == ".":
                    alt = "-"
                else:
                    alt = info[2]
                refDict[refcontig][position-1] = alt
    # write output fasta file
    with open(iso+".fasta","w") as out:
        for contig,seq in refDict.items():
            seqstr = "".join(seq)
            chunks = [seqstr[i:i+60] for i in range(0,len(seqstr),60)]
            out.write(">"+iso+"_"+contig+"\n")
            for chunk in chunks:
                out.write(chunk+"\n")

for fasta in fastas:
    runMummer(fasta,reffile)
