#!/usr/bin/env python3

import sys
import os
sys.path.append("../helper")
from myio import get_strain_lst
from gff import read_gff

target = sys.argv[1]
strain_lst = get_strain_lst(target, full = True)
direc = "/data/mitsuki/data/denovo/{}".format(target)
for strain in strain_lst:
    dnafp = "{}/dnaseq/{}.dnaseq".format(direc, strain)
    fnafp = "{}/annotation/refseq/fna/{}.fna".format(direc, strain)
    faafp = "{}/annotation/refseq/faa/{}.faa".format(direc, strain)
    gfffp = "{}/annotation/refseq/gff/{}.gff".format(direc, strain)

    for fp in (dnafp, fnafp, faafp):
        if not(os.path.isfile(fp)):
            print("ERROR: {} not found".format(fp), file = sys.stderr)
            exit(1)

    gff_df = read_gff(gfffp, ["orf_id"])
    if not("orf_id" in gff_df.columns):
        print("ERROR: orf_id not found in {}".format(gfffp), file = sys.stderr)
        exit(2)
        
exit(0)
