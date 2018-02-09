#!/usr/bin/env python3

"""
prepare input directory for sonicparanod will run on
"""

import sys
import os
from Bio import SeqIO

sys.path.append("../helper")
from myio import get_strain_lst

target=sys.argv[1]
annotationType=sys.argv[2]
strain_lst = get_strain_lst(target)

baseDirec="/data/mitsuki/data/denovo/{}/annotation/{}".format(target, annotationType)
os.makedirs("{}/sonic/in".format(baseDirec), exist_ok=True)
os.makedirs("{}/sonic/out".format(baseDirec), exist_ok=True)

print("START: prepare {} input faa".format(len(strain_lst)))
for strain in strain_lst:
    inFilepath="{}/faa/{}.faa".format(baseDirec, strain)
    outFilepath="{}/sonic/in/{}".format(baseDirec, strain)  #sonicparanoid can not handle extention
    with open(outFilepath, "w") as f:
        for rec in SeqIO.parse(inFilepath, "fasta"):
            rec.description="" #sonicparanoid can not handle description column
            SeqIO.write(rec, f, "fasta")
print("DONE: prepare {}/sonic".format(baseDirec))
    
