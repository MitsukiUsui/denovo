#!/usr/bin/env python3

import sys
import pandas as pd
from Bio import SeqIO

sys.path.append("../helper")
from myio import get_strain_lst
from gff import read_gff
from interval import Interval, interval_sum, interval_justsum, interval_not, interval_and


def calc_stat(target, strain):
    seqFilepath = "/data/mitsuki/data/denovo/{}/dnaseq/{}.dnaseq".format(target, strain)
    rec_lst = []
    for rec in SeqIO.parse(seqFilepath, "fasta") :
        rec_lst.append(rec)

    hitFilepath = "../blastn/result/{}/{}.csv".format(target, strain)
    hit_df = pd.read_csv(hitFilepath)

    gffFilepath = "/data/mitsuki/data/denovo/{}/annotation/refseq/gff/{}.gff".format(target, strain)
    gff_df = read_gff(gffFilepath, ["pseudo"])
    if "pseudo" in gff_df.columns:
        gff_df["pseudo_b"] = gff_df["pseudo"].notnull()
    else:
        gff_df["pseudo_b"] = False
        
    length = 0
    genicSum = 0
    pseudoSum = 0
    interSum = 0
    hitSum = 0
    hitJustSum = 0
    genicHitSum = 0
    pseudoHitSum= 0
    interHitSum = 0
    
    for rec in rec_lst:
        seqname = rec.id

        # define genic_lst, psedo_lst, inter_lst according to gff_df
        genic_lst = []
        pseudo_lst = []
        for _, row in gff_df[gff_df["seqname"] == seqname].iterrows():
            if row["pseudo_b"]:
                pseudo_lst.append( Interval(row["start"] -1, row["end"]) )
            else:
                genic_lst.append( Interval(row["start"] -1, row["end"]) )
        inter_lst =  interval_not(genic_lst + pseudo_lst, len(rec))

        # define hit_lst according to hit_df
        hit_lst = []
        for _, row in hit_df[hit_df["sseqid"] == seqname].iterrows():
            if row["hit_strand"] == 1:
                start = row["sstart"] - 1
                end = row["send"]
            else:
                start = row["send"] - 1
                end = row["sstart"]
            hit_lst.append(Interval(start, end))

        length += len(rec)
        genicSum += interval_sum(genic_lst)
        pseudoSum += interval_sum(pseudo_lst)
        interSum += interval_sum(inter_lst)
        hitSum += interval_sum(hit_lst)
        hitJustSum += interval_justsum(hit_lst)
        genicHitSum += interval_and(genic_lst, hit_lst)
        pseudoHitSum += interval_and(pseudo_lst, hit_lst)
        interHitSum += interval_and(inter_lst, hit_lst)
        
    ret = {
        "length": length,
        "genic_sum" : genicSum,
        "pseudo_sum": pseudoSum,
        "inter_sum": interSum,
        "hit_sum": hitSum,
        "hit_justsum": hitJustSum,
        "genic_hit_sum": genicHitSum,
        "pseudo_hit_sum": pseudoHitSum,
        "inter_hit_sum": interHitSum
    }
    return ret


def main(target, outFilepath):
    strain_lst = get_strain_lst(target)
    print("START: collect blastn stat for {} strains".format(len(strain_lst)), flush = True)
    
    dct_lst = []
    for strain in strain_lst:
        dct = {"strain": strain}
        dct.update(calc_stat(target, strain))
        dct_lst.append(dct)
    stat_df = pd.DataFrame(dct_lst)
    stat_df = stat_df[["strain", "length", "genic_sum", "pseudo_sum", "inter_sum", 
                            "hit_sum", "hit_justsum",  "genic_hit_sum", "pseudo_hit_sum", "inter_hit_sum"]]
    stat_df.to_csv(outFilepath, index=False)
    print("DONE: output {}".format(outFilepath), flush = True)

if __name__=="__main__":
    target = sys.argv[1]
    outFilepath = "./out/blastn/{}.csv".format(target)
    main(target, outFilepath)
