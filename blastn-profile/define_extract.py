#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

sys.path.append("../helper")
from gff import read_gff

def get_query_df(region_df):
    """
    :ret qof_id, qo(start/end)_(dna/pro), qstrand
    """

    ret_df = region_df.rename(columns={"rqstart": "qostart_dna", "rqend": "qoend_dna"})
    ret_df["qostart_pro"] = np.ceil(ret_df["qostart_dna"] / 3).astype(int)
    ret_df["qoend_pro"] = np.floor(ret_df["qoend_dna"] / 3).astype(int)
    ret_df = ret_df[["region_id", "qorf_id", "qostart_dna", "qoend_dna", "qostart_pro", "qoend_pro", "qstrand"]]
    return ret_df
    
    
def get_sbjct_df(region_df, gff_df):
    """
    :ret sof_id, so(start/end)_(dna/pro), sstrand
    """
    
    gff_df = gff_df.rename(columns={"orf_id":"sorf_id", "start": "cstart", "end": "cend", "strand": "sstrand"})
    gff_df["cstart"] = gff_df["cstart"] - 1 # c for cds
    gff_df["sstrand"] = [1 if (strand == "+") else -1 for strand in gff_df["sstrand"]]
    join_df = pd.merge(region_df, gff_df[["sorf_id", "cstart", "cend", "sstrand"]], on="sorf_id", how="left")

    sostartDna_lst = []
    soendDna_lst = []
    for _, row in join_df.iterrows():
        if row["sstrand"] == 1:
            sostartDna_lst.append(row["ostart"] - row["cstart"])
            soendDna_lst.append(row["oend"] - row["cstart"])
        elif row["sstrand"] == -1:
            sostartDna_lst.append(row["cend"] - row["oend"])
            soendDna_lst.append(row["cend"] - row["ostart"])
    join_df["sostart_dna"] = sostartDna_lst
    join_df["soend_dna"] = soendDna_lst
    join_df["sostart_pro"] = np.ceil(join_df["sostart_dna"] / 3).astype(int)
    join_df["soend_pro"] = np.floor(join_df["soend_dna"] / 3).astype(int)

    ret_df = join_df[["region_id", "sorf_id", "sostart_dna", "soend_dna", "sostart_pro", "soend_pro", "sstrand"]]
    return ret_df

def main(regionFilepath, gffFilepath, outFilepath):
    region_df = pd.read_csv(regionFilepath)
    region_df = region_df[region_df["category"] == "genic"] #filter only "genic"
    
    query_df = get_query_df(region_df)
    
    gff_df = read_gff(gffFilepath, ["orf_id"])
    sbjct_df = get_sbjct_df(region_df, gff_df)

    out_df = pd.merge(query_df, sbjct_df, on="region_id")
    out_df = pd.merge(out_df, region_df[["region_id", "olength"]], on="region_id")
    assert out_df.shape[0] == region_df.shape[0]

    out_df.to_csv(outFilepath, index=False)
    print("DONE: output {}".format(outFilepath))

if __name__ == "__main__":
    target=sys.argv[1]
    strain=sys.argv[2]
    
    direc = "./out/{}/{}".format(target, strain)
    regionFilepath = "{}/region.csv".format(direc)
    gffFilepath = "/data/mitsuki/data/denovo/{}/annotation/refseq/gff/{}.gff".format(target, strain)
    outFilepath = "{}/extract.csv".format(direc)
    main(regionFilepath, gffFilepath, outFilepath)
