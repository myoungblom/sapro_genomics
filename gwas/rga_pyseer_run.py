#!/usr/bin/env/python3

import sys
import subprocess

tests = ["human","animal","food","builtenv","naturalenv"]

vcffile = "masked_ncbiAlns_rgas_snps.vcf"
for t in tests:
    print(t)
    phenofile = "rga_phenotypes/sapro_"+t+"_pyseer_traits.txt"
    patterns = "allRGA_"+t+"_patterns.txt"
    outfile = "sapro_allRGA_"+t+"_MixedModel.txt"
    with open(outfile,"w") as output:
        subprocess.call(["pyseer","--phenotypes",phenofile,"--vcf",vcffile,"--lmm","--similarity","sapro_RGA_distances.tsv",\
                    "--output-patterns",patterns],stdout=output)
