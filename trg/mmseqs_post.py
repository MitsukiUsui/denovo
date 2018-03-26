#!/usr/bin/env python3

"""
append family, tax_id, qlength information to .m8
"""

import sys
import pandas as pd
from Bio import SeqIO

sys.path.append("../helper")
from myio import get_cluster_df, get_strain_lst
from TaxDbController import TaxDbController

def read_m8(fp):
    columns_lst=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart","send", "evalue", "bitscore"]
    hit_df=pd.read_csv(fp, delimiter='\t', header=None, names=columns_lst)
    return hit_df

def get_taxid_lst(accessionVersion_lst):
    tdc = TaxDbController("/data/mitsuki/data/refseq/accession2taxid/prot.accession2taxid.sq3")

    #  create accession2taxid lookup dictionaly first
    accessionVersion_set=set(accessionVersion_lst)
    accession2taxid=dict()
    print("START: assign taxid to {} unique accession".format(len(accessionVersion_set)), flush=True)
    batch=len(accessionVersion_set)*0.01 #1% of total lines
    border=batch
    for _,accessionVersion in enumerate(accessionVersion_set):
        if _>=border:
            border+=batch
            print(".", end="", flush=True)
        accession2taxid[accessionVersion]=tdc.accession2taxid(accessionVersion)
    print()

    taxid_lst=[accession2taxid[accessionVersion] for accessionVersion in accessionVersion_lst]
    return taxid_lst

def get_qlength_lst(orfid_lst, queryFilepath):
    id2length = {}
    for rec in SeqIO.parse(queryFilepath, "fasta"):
        id2length[rec.id] = len(rec.seq)
    qlength_lst = [id2length[orfid] for orfid in orfid_lst]
    return qlength_lst

def main(target, queryFilepath, resultFilepath, outFilepath):
    result_df = read_m8(resultFilepath)
    result_df = result_df.rename(columns={'qseqid': 'orf_id', 'sseqid': "accession_version"})
    print("DONE: load {} hit result {}".format(result_df.shape[0], resultFilepath), flush=True)

    # create family_df ("family", "orf_id")
    strain_lst = get_strain_lst(target)
    cluster_df = get_cluster_df(target)
    dct_lst = []
    for _, row in cluster_df.iterrows():
        family = row["family"]
        for orfids in row[strain_lst].dropna():
            for orfid in orfids.split(","):
                dct = {"family": family, "orf_id": orfid}
                dct_lst.append(dct)
    family_df = pd.DataFrame(dct_lst)

    join_df = pd.merge(result_df, family_df, on="orf_id", how="left")
    print("DONE: merge family", flush=True)
    join_df["tax_id"] = get_taxid_lst(join_df["accession_version"])
    print("DONE: merge taxonomy", flush=True)
    join_df["qlength"] = get_qlength_lst(join_df["orf_id"], queryFilepath)
    print("DONE: merge qlength", flush=True)

    col_lst = ["family", "orf_id", "accession_version", "tax_id",
                "pident", "length", "mismatch", "gapopen", "qstart", "qend", "qlength", "sstart", "send", "evalue", "bitscore"]
    join_df = join_df[col_lst]
    join_df.to_csv(outFilepath, index=False)
    print("DONE: output {}".format(outFilepath), flush=True)

if __name__ == "__main__":
    target = sys.argv[1]
    queryFilepath = sys.argv[2]
    resultFilepath = sys.argv[3]
    outFilepath = sys.argv[4]
    main(target, queryFilepath, resultFilepath, outFilepath)

