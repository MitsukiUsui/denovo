#!/usr/bin/env python3

import math
import sys
import os
import numpy as np
import pandas as pd

sys.path.append("../helper")
from gff import read_gff

def get_overlap_df(gff_df, hit_df):
    def get_overlap_dctdct(gff_df, hit_df, chrName):

        # 1. create ordered position of list together with CDS and blastn hits(regions).
        #    format: (pos, isStart, type, id)
        #    * type: 0=CDS, 1=region
        
        filtered_df=hit_df[hit_df["sseqid"]==chrName]
        pos_lst=[]
        for _, row in gff_df.iterrows():
            if row["seqname"]==chrName and row["feature"]=="CDS":
                pos_lst.append((row["start"] - 1, True,  0, row["orf_id"]))
                pos_lst.append((row["end"], False, 0, row["orf_id"]))
        for _, row in hit_df.iterrows():
            if row["sseqid"]==chrName:
                stt = min(row["sstart"], row["send"]) - 1
                end = max(row["sstart"], row["send"])
                pos_lst.append( (stt, True, 1, row["hit_id"]) )
                pos_lst.append( (end, False, 1, row["hit_id"]) )
        pos_lst=sorted(pos_lst, key=lambda x: (x[0], x[1], x[2])) # sort by first element. need to be deterministic...?

        # 2. create overlap_dctdct
        #  key: (orf_id, hit_id), value: (ofirst, olast)
        
        overlap_dctdct={}
        cds_lst, region_lst=[], []  # remember in process cds or region
        for pos in pos_lst:
            if pos[2]==0:#if CDS
                for region in region_lst:
                    key = (pos[3],  region)
                    if pos[1]: #if start initialize new dct
                        overlap_dctdct[key]={}
                        overlap_dctdct[key]["ostart"] = pos[0]
                    else: # if end 
                        overlap_dctdct[key]["oend"] = pos[0]
                if pos[1]:
                    cds_lst.append(pos[3])
                else:
                    cds_lst.remove(pos[3])
            elif pos[2]==1:#if region
                for cds in cds_lst:
                    key = (cds, pos[3])
                    if pos[1]:
                        overlap_dctdct[key]={}
                        overlap_dctdct[key]["ostart"]=pos[0]
                    else:
                        overlap_dctdct[key]["oend"]=pos[0]
                if pos[1]:#if start
                    region_lst.append(pos[3])
                else:
                    region_lst.remove(pos[3])                        
        return overlap_dctdct
    
    chrName_lst=list(set(hit_df["sseqid"]))
    dct_lst=[]
    for chrName in chrName_lst:
        print("START: process {}".format(chrName))
        overlap_dctdct=get_overlap_dctdct(gff_df, hit_df, chrName)
        
        #update dct_lst according to overlap_dctdct
        for k,v in overlap_dctdct.items():
            dct={}
            dct["sorf_id"], dct["hit_id"] = k[0], k[1]
            dct["ostart"], dct["oend"] = v["ostart"], v["oend"]
            dct["olength"]=v["oend"] - v["ostart"]
            dct["chr_name"]=chrName
            dct_lst.append(dct)
        print("\tfound {} overlaps".format(len(overlap_dctdct)))
        
    overlap_df=pd.DataFrame(dct_lst)
    overlap_df["overlap_id"] = overlap_df.index
    overlap_df=overlap_df[["overlap_id", "sorf_id", "hit_id", "chr_name", "ostart", "oend", "olength"]]
    return overlap_df

def add_sbjct_pos(overlap_df, gff_df):
    # format gff_df & left join to overlap_df
    gff_df = gff_df.rename(columns={"orf_id":"sorf_id"})
    gff_df["cstart"] = gff_df["start"] - 1 # c for cds
    gff_df["cend"] = gff_df["end"]
    gff_df["sstrand"] = [1 if (strand == "+") else -1 for strand in gff_df["strand"]]
    gff_df["sfamily"] = gff_df["family"]
    gff_df=gff_df[["sorf_id", "sstrand", "sfamily", "cstart", "cend"]]
    ret_df = pd.merge(overlap_df, gff_df, on = "sorf_id", how="left")
    assert ret_df.shape[0] == overlap_df.shape[0]

    # calculate sostart & soend based on cfirst and clast
    dct_lst=[] 
    for _, row in ret_df.iterrows():
        dct={}
        if row["sstrand"] == 1:
            dct["sostart_dna"] = row["ostart"] - row["cstart"]
            dct["soend_dna"] = row["oend"] - row["cstart"]
        elif row["sstrand"] == -1:
            dct["sostart_dna"] = row["cend"] - row["oend"]
            dct["soend_dna"] = row["cend"] - row["ostart"]
        dct["sostart_pro"]=math.ceil(dct["sostart_dna"] / 3)
        dct["soend_pro"]=math.floor(dct["soend_dna"] / 3)
        dct_lst.append(dct)
    tmp_df=pd.DataFrame(dct_lst, index=ret_df.index)
    ret_df=pd.concat([ret_df, tmp_df], axis=1)
    return ret_df

def add_query_pos(overlap_df, hit_df):
    # format hit_df & left join to overlap_df
    hit_df["qstrand"] = hit_df["hit_strand"]
    hit_df["qorf_id"] = hit_df["qseqid"]
    hit_df=hit_df[["hit_id", "qstart", "qend", "sstart", "send", "qstrand", "qorf_id", "qfamily", "qstrain", "sstrain"]]
    ret_df=pd.merge(overlap_df, hit_df, on="hit_id")
    assert ret_df.shape[0]==overlap_df.shape[0]
    
    # calculate qostart & qoend
    dct_lst=[]
    for _,row in ret_df.iterrows():
        #calc queryOvelapStartPercentage(qosp) & queryOverlapEndPercentage(qoep)
        if row["qstrand"] == 1:
            length = row["send"] - row["sstart"] + 1
            assert length > 0
            qosp = ( row["ostart"] - (row["sstart"] - 1) ) / length  #-1 to convert to 0 origin half-open interval
            qoep = ( row["oend"] - (row["sstart"] - 1) ) / length
        elif row["qstrand"] == -1:
            length = row["sstart"] - row["send"] + 1
            assert length > 0
            qosp = ( row["sstart"] - row["oend"]) / length  #we can treat row["sstart"] as the latter of 0 origin half-open interval
            qoep = ( row["sstart"] - row["ostart"]) / length
        else:
            print("ERROR: undefined strand {}".format(row["hit_strand"]))
            sys.exit()
        
        # calc qostart & qoend based on qosp & qoep each
        length = row["qend"] - row["qstart"] + 1
        qostart_dna = math.ceil(row["qstart"] - 1 + length * qosp)
        qoend_dna = math.floor(row["qstart"] - 1 + length * qoep)
        qostart_pro = math.ceil(qostart_dna / 3)
        qoend_pro = math.floor(qoend_dna / 3)
        
        dct = {"qosp": qosp, "qostart_dna": qostart_dna, "qostart_pro": qostart_pro, 
                   "qoep": qoep, "qoend_dna": qoend_dna, "qoend_pro": qoend_pro}
        dct_lst.append(dct)
    tmp_df=pd.DataFrame(dct_lst, index=ret_df.index)
    ret_df=pd.concat([ret_df, tmp_df], axis=1)
    return ret_df

def main(strain, hitFilepath, geneFilepath, overlapFilepath):
    hit_df=pd.read_csv(hitFilepath)
    gff_df=read_gff(gffFilepath, ["orf_id", "family"])
   
    ovr_df=get_overlap_df(gff_df, hit_df)
    ovr_df=add_sbjct_pos(ovr_df, gff_df)
    ovr_df=add_query_pos(ovr_df, hit_df)
    
    column_lst=["overlap_id", "hit_id", "ostart", "oend", "olength", "chr_name",
                "qstrain", "sstrain", "qfamily", "sfamily", "qstrand", "sstrand", "qorf_id", "sorf_id",
                "qostart_dna", "qostart_pro", "qoend_dna", "qoend_pro",
                "sostart_dna", "sostart_pro", "soend_dna", "soend_pro",
                "qstart", "qend", "sstart", "send", "cstart", "cend", "qosp", "qoep"]
    ovr_df = ovr_df[column_lst]
    ovr_df = ovr_df.set_index("overlap_id")
    ovr_df.to_csv(overlapFilepath)
    print("DONE: {} overlaps in {}".format(ovr_df.shape[0], overlapFilepath))

if __name__=="__main__":
    target=sys.argv[1]
    strain=sys.argv[2]
    annotationType="refseq"
    
    hitFilepath="../blastn/result/{}/{}.csv".format(target, strain)
    gffFilepath="/data/mitsuki/data/denovo/{}/annotation/{}/gff/{}.gff".format(target, annotationType, strain)
    outDirec="./out/{}".format(target)
    os.makedirs(outDirec, exist_ok=True)
    overlapFilepath="{}/{}_ovr.csv".format(outDirec, strain)
    main(strain,hitFilepath, gffFilepath, overlapFilepath)
