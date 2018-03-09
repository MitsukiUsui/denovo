#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
from Bio import SeqIO

sys.path.append("../helper")
from gff import read_gff
from myutil import myinterval

def format_out_df(out_df, hit_df):
    """
    extract query information from hit_df
    append rqstart, rqend, rstrand, qorf_id
    """
    
    #calculate rqstart, rqend
    hit_df = hit_df.rename(columns={"qseqid" : "qorf_id", "hit_strand" : "qstrand"})
    join_df = pd.merge(out_df, hit_df[["hit_id", "qstart", "qend", "sstart", "send", "qorf_id", "qstrand"]], on="hit_id", how="left")
    qosp_lst = []
    qoep_lst = []
    for _, row in join_df.iterrows():
        if row["qstrand"] == 1:
            length = row["send"] - row["sstart"] + 1
            assert length > 0
            qosp = ( row["ostart"] - (row["sstart"] - 1) ) / length  #-1 to convert to 0 origin half-open interval
            qoep = ( row["oend"] - (row["sstart"] - 1) ) / length
        elif row["qstrand"] == -1:
            length = row["sstart"] - row["send"] + 1
            assert length > 0
            qosp = ( row["sstart"] - row["oend"]) / length  #we can treat row["sstart"] (1 origin closed end) as the latter of 0 origin open end
            qoep = ( row["sstart"] - row["ostart"]) / length
        else:
            print("ERROR: undefined strand {}".format(row["hit_strand"]), file=sys.stderr)
            sys.exit()
        qosp_lst.append(qosp)
        qoep_lst.append(qoep)
    join_df["qosp"] = qosp_lst
    join_df["qoep"] = qoep_lst
    join_df["rqstart"] = np.ceil(join_df["qstart"]-1 + (join_df["qend"]-join_df["qstart"]+1) * join_df["qosp"]).astype(int)
    join_df["rqend"] = np.floor(join_df["qstart"]-1 + (join_df["qend"]-join_df["qstart"]+1) * join_df["qoep"]).astype(int)

    ret_df = join_df[list(out_df.columns) + ["rqstart", "rqend", "qstrand", "qorf_id"]]
    return ret_df

def get_seq2length(dnaFilepath):
    seq2length = {}
    for rec in SeqIO.parse(dnaFilepath, "fasta"):
        seq2length[rec.id] = len(rec.seq)
    return seq2length

def main(hitFilepath, gffFilepath, dnaFilepath, outFilepath):
    hit_df = pd.read_csv(hitFilepath)
    
    gff_df = read_gff(gffFilepath, ["orf_id","pseudo"])
    if "pseudo" in gff_df.columns:
        gff_df["pseudo_b"] = gff_df["pseudo"].notnull()
    else:
        gff_df["pseudo_b"] = False

    dct_lst = []
    seq2length = get_seq2length(dnaFilepath)
    for seqname, seqlength in seq2length.items():
        #define list of Interval for genic, pseudogene, intergenic regions respectively
        genic_lst = []
        pseudo_lst = []
        for _, row in gff_df[gff_df["seqname"] == seqname].iterrows():
            if row["pseudo_b"]:
                pseudo_lst.append(myinterval.Interval(row["start"]-1, row["end"], row["orf_id"]))
            else:
                genic_lst.append(myinterval.Interval(row["start"]-1, row["end"], row["orf_id"]) )
        inter_lst =  myinterval.complement(genic_lst+pseudo_lst, 0,  seqlength, setid=True)
        
        # define list of Interval for hit
        hit_lst = []
        for _, row in hit_df[hit_df["sseqid"] == seqname].iterrows():
            if row["hit_strand"] == 1:
                hit_lst.append(myinterval.Interval(row["sstart"]-1, row["send"], row["hit_id"]))
            else:
                hit_lst.append(myinterval.Interval(row["send"]-1, row["sstart"], row["hit_id"]))
            
        # calcurate intersection 
        hitGenic_lst = myinterval.intersection(hit_lst, genic_lst)
        hitPseudo_lst = myinterval.intersection(hit_lst, pseudo_lst)
        hitInter_lst = myinterval.intersection(hit_lst, inter_lst)
        
        #TODO: add seqname to dict
        for interval in hitGenic_lst:
            dct = {"hit_id": interval.id1,
                   "sorf_id": interval.id2,
                   "ostart": interval.start,
                   "oend": interval.end,
                   "seqname": seqname,
                   "category": "genic"}
            dct_lst.append(dct)

        for interval in hitPseudo_lst:
            dct = {"hit_id": interval.id1,
                   "sorf_id": interval.id2,
                   "ostart" : interval.start,
                   "oend": interval.end,
                   "seqname": seqname,
                   "category": "pseudo"}
            dct_lst.append(dct)

        for interval in hitInter_lst:
            dct = {"hit_id": interval.id1,
                   "ostart": interval.start,
                   "oend": interval.end,
                   "seqname": seqname,
                   "category": "inter"}
            dct_lst.append(dct) 
        
    out_df = pd.DataFrame(dct_lst)
    out_df["olength"] = out_df["oend"] - out_df["ostart"]
    out_df = format_out_df(out_df, hit_df)
    out_df = out_df[["hit_id", "category", "seqname", "sorf_id", "ostart", "oend", "olength", "qorf_id", "rqstart", "rqend", "qstrand"]]
    out_df = out_df.sort_values(by=["hit_id", "ostart"])
    out_df = out_df.reset_index(drop=True)
    out_df.index.name = "region_id"
    out_df.to_csv(outFilepath)
    print("DONE: output {}".format(outFilepath))

if __name__=="__main__":
    target = sys.argv[1]
    strain=sys.argv[2]

    hitFilepath = "../blastn/result/{}/{}.csv".format(target, strain)
    gffFilepath = "/data/mitsuki/data/denovo/{}/annotation/refseq/gff/{}.gff".format(target, strain)
    dnaFilepath = "/data/mitsuki/data/denovo/{}/dnaseq/{}.dnaseq".format(target, strain)
    outFilepath = "./out/{}/{}/region.csv".format(target, strain)
    main(hitFilepath, gffFilepath, dnaFilepath, outFilepath)
