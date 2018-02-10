#!/usr/bin/env python3
import sys

sys.path.append("../helper")
from myio import get_strain_lst

target=sys.argv[1]
baseDirec="/data/mitsuki/data/denovo/{}".format(target)
strain_lst = get_strain_lst(target)
for strain in strain_lst:
    seqFilepath="{}/dnaseq/{}.dnaseq".format(baseDirec, strain)
    gffFilepath="{}/annotation/prodigal/gff/{}.gff".format(baseDirec, strain)
    fnaFilepath="{}/annotation/prodigal/fna/{}.fna".format(baseDirec, strain)
    faaFilepath="{}/annotation/prodigal/faa/{}.faa".format(baseDirec, strain)
    supFilepath="{}/annotation/prodigal/sup/{}.sup".format(baseDirec, strain)
    print("{},{},{},{},{}".format(seqFilepath, gffFilepath, fnaFilepath, faaFilepath, supFilepath))

