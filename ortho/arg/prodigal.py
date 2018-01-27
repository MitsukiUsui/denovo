#!/usr/bin/env python3
import sys

target=sys.argv[1]
baseDirec="/data/mitsuki/data/denovo/{}".format(target)
strainFilepath="../data/{}/strain.lst".format(target)
strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
for strain in strain_lst:
    seqFilepath="{}/dnaseq/{}.dnaseq".format(baseDirec, strain)
    gffFilepath="{}/annotation/prodigal/gff/{}.gff".format(baseDirec, strain)
    fnaFilepath="{}/annotation/prodigal/fna/{}.fna".format(baseDirec, strain)
    faaFilepath="{}/annotation/prodigal/faa/{}.faa".format(baseDirec, strain)
    print("{},{},{},{}".format(seqFilepath, gffFilepath, fnaFilepath, faaFilepath))


