#!/usr/bin/env python3

"""
prepare input directory for sonicparanod will run on
"""

import sys
import os
from Bio import SeqIO

target=sys.argv[1]
annotationType=sys.argv[2]

strainFilepath="../data/{}/strain.lst".format(target)
strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
baseDirec="/data/mitsuki/data/denovo/{}/annotation/{}".format(target, annotationType)
sonicDirec="{}/sonic".format(baseDirec)
os.makedirs("{}/in".format(sonicDirec), exist_ok=True)
os.makedirs("{}/out".format(sonicDirec), exist_ok=True)

print("START: prepare {} input faa".format(len(strain_lst)))
for strain in strain_lst:
    inFilepath="{}/faa/{}.faa".format(baseDirec, strain)
    outFilepath="{}/in/{}".format(sonicDirec, strain)  #sonicparanoid can not handle extention
    with open(outFilepath, "w") as f:
        for rec in SeqIO.parse(inFilepath, "fasta"):
            rec.description="" #sonicparanoid can not handle description column
            SeqIO.write(rec, f, "fasta")
print("DONE: prepare {}".format(sonicDirec))
    
