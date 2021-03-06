#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO
from collections import defaultdict

sys.path.append("../helper")
from gff import read_gff
from myio import *

def output_family2rec(family2rec, outDirec, ext):
    for family, rec_lst in family2rec.items():
        outFilepath="{}/{}.{}".format(outDirec, family, ext)
        with open(outFilepath, "w") as f:
            for rec in rec_lst:
                SeqIO.write(rec, f, "fasta")

def main(target, strain_lst):
    annotType="refseq"
    direc="/data/mitsuki/data/denovo/{}/annotation/{}".format(target, annotType)
    
    print("START: parse {} * 2 FASTA files".format(len(strain_lst)))
    family2fna=defaultdict(list)
    family2faa=defaultdict(list)
    for strain in strain_lst:
        gffFilepath="{}/gff/{}.gff".format(direc, strain)
        fnaFilepath="{}/fna/{}.fna".format(direc, strain)
        faaFilepath="{}/faa/{}.faa".format(direc, strain)

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
    outDirec="{}/family/fna".format(direc)
    os.makedirs(outDirec, exist_ok=True)
    output_family2rec(family2fna, outDirec, "fna")
    print("\tDONE: output {} family in {}".format(len(family2fna), outDirec))

    outDirec="{}/family/faa".format(direc)
    os.makedirs(outDirec, exist_ok=True)
    output_family2rec(family2faa, outDirec, "faa")
    print("\tDONE: output {} family in {}".format(len(family2faa), outDirec))

if __name__=="__main__":
    target=sys.argv[1]
    strain_lst = get_strain_lst(target)
    main(target, strain_lst)
