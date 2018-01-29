#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO
from collections import defaultdict

sys.path.append("../helper")
from gff import read_gff

def output_family2rec(family2rec, outDirec, ext):
    for family, rec_lst in family2rec.items():
        outFilepath="{}/{}.{}".format(outDirec, family, ext)
        with open(outFilepath, "w") as f:
            for rec in rec_lst:
                SeqIO.write(rec, f, "fasta")

def main(target, strain_lst):
    prodigalDirec="/data/mitsuki/data/denovo/{}/annotation/prodigal".format(target)
    
    print("START: parse {} * 2 FASTA files".format(len(strain_lst)))
    family2fna=defaultdict(list)
    family2faa=defaultdict(list)
    for strain in strain_lst:
        gffFilepath="{}/gff/{}.gff".format(prodigalDirec, strain)
        fnaFilepath="{}/fna/{}.fna".format(prodigalDirec, strain)
        faaFilepath="{}/faa/{}.faa".format(prodigalDirec, strain)

        gff_df=read_gff(gffFilepath, ["orf_id","family"])
        id2family={}
        for _, row in gff_df.iterrows():
            id2family[row["orf_id"]]=row["family"]

        for rec in SeqIO.parse(fnaFilepath, "fasta"):
            family=id2family[rec.id]
            family2fna[family].append(rec)

        for rec in SeqIO.parse(faaFilepath, "fasta"):
            family=id2family[rec.id]
            family2faa[family].append(rec)

    print("START: output fna & faa for every family")
    outDirec="{}/family/fna".format(prodigalDirec)
    os.makedirs(outDirec, exist_ok=True)
    output_family2rec(family2fna, outDirec, "fna")
    print("\tDONE: output {} family in {}".format(len(family2fna), outDirec))

    outDirec="{}/family/faa".format(prodigalDirec)
    os.makedirs(outDirec, exist_ok=True)
    output_family2rec(family2faa, outDirec, "faa")
    print("\tDONE: output {} family in {}".format(len(family2faa), outDirec))

if __name__=="__main__":
    target=sys.argv[1]
    strainFilepath="../data/{}/strain.lst".format(target)
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    main(target, strain_lst)