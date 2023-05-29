#!/usr/bin/env/python3

import argparse
import sys
import os
import getopt
import egglib
print("Egglib version: "+egglib.__version__)
import numpy as np
from Bio import SeqIO
from collections import defaultdict

####
# This script reads an alignment and calculates diversity and selection statistics
# in a sliding window over a subsample of the alignment, repeated 100x.
####

def get_args():
    """
    Handle the command line arguments.
    """
    parser = argparse.ArgumentParser(description="This script reads an alignment and calculates diversity and selection statistics\
    in a sliding window over a subsample of the alignment, repeated 100x.")
    parser.add_argument("alnfile",help="alignment file in fasta format")
    parser.add_argument("-w",default=20000,type=int,help="window size. default = 20,000 bp")
    parser.add_argument("-s",default=5000,type=int,help="step size. default= 5,000 bp")
    parser.add_argument("-n",default=50,type=int,help="sub-sample size. default = 50")
    parser.add_argument("-x",default=100,type=int,help="number of iterations. default = 100")
    return parser.parse_args()

args = get_args()

def readaln(alnfile):
    """
    Read alignment file into a list.
    """
    seqs = []
    for seq in SeqIO.parse(alnfile,"fasta"):
        seqs.append(seq)
    return seqs

def subsample(seqlist,n,outfile):
    """
    Sub samples list of sequences to size n.
    """
    indices = list(np.random.choice(range(len(seqlist)), size=n, replace=False))
    sample = [seqlist[i] for i in indices]
    SeqIO.write(sample,outfile,"fasta")
    outfile.close()

# read full alignment, write 100 subsampled alignments
print("Reading alignment ...")
all_alns = readaln(args.alnfile)
aln_name = args.alnfile.split("_")[0]
print("Subsampling alignment ...")
for i in range(1,101):
    outfile = open(aln_name+"_subsample_"+str(i)+".fasta","w")
    subsample(all_alns,args.n,outfile)
    outfile.close()

def calc_stats(aln):
    """
    Calculate statistics for alignment.
    """
    statDict = {}
    cs = egglib.stats.ComputeStats(multi_hits=True)
    cs.add_stats('thetaW', 'Pi', 'D','lseff')
    polyDict = cs.process_align(aln,max_missing=0.15)
    try:
        statDict['theta'] = polyDict['thetaW']/polyDict['lseff']
        statDict['pi'] = polyDict['Pi']/polyDict['lseff']
        statDict['tajimasD'] = polyDict['D']
    except TypeError:
        statDict['theta'] = "NA"
        statDict['pi'] = "NA"
        statDict['tajimasD'] = "NA"
    return statDict


def calcSubsampleStats(width,step):
    """
    Calculate & store stats for all subsamples.
    """
    subsample_dict = defaultdict(list)
    for i in range(1,101):
        start = 0
        stop = width
        print("subsample "+str(i))
        infile = aln_name+"_subsample_"+str(i)+".fasta"
        align = egglib.io.from_fasta(infile,alphabet=egglib.alphabets.DNA)
        for window in align.slider(width,step):
            stats = calc_stats(window)
            location = int((start+stop)/2)
            subsample_dict[location].append([stats['theta'],stats['pi'],stats['tajimasD']])
            start += step
            stop += step
    return subsample_dict

print("Calculating sliding window stats ...")
subDict = calcSubsampleStats(args.w,args.s)

print("Writing output ...")
outfile = open(aln_name+'windowStats_w'+str(args.w)+"_s"+str(args.s)+'.txt', 'w')
outfile.write("Position\tTheta\tPi\tTajimasD\n")
for position,stats in subDict.items():
    for sub in stats:
        outfile.write("\t".join([str(position)]+[str(x) for x in sub])+"\n")
