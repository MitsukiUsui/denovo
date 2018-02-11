#!/usr/bin/env python3

import sys
from Bio import SeqIO
import pandas as pd
#import numpy as np
#from Bio.Seq import Seq

def extract_sequence(orf2record, orfId, start, end):
    return orf2record[orfId].seq[start:end]

def main(target, ovrFilepath, seqFilepath, orf2fna, orf2faa):
    overlap_df=pd.read_csv(ovrFilepath, index_col=0)
    msk = (overlap_df["olength"] > 10) & (~overlap_df["sfamily"].isnull()) # extract only when overlap length > 10, and sbjct is not pseudogene
    overlap_df = overlap_df[msk]
    
    print("START: extract sequence for {} / {} overlaps".format(sum(msk), len(msk)), flush = True)
    with open(seqFilepath, 'w') as f:
        for key,row in overlap_df.iterrows():
            if (row["sorf_id"] not in orf2fna) or (row["sorf_id"] not in orf2faa):
                print("ERROR: {} does not have sequence information".format(row["sorf_id"]), file = sys.stderr)
            elif (row["qorf_id"] not in orf2fna) or (row["qorf_id"] not in orf2faa):
                print("ERROR: {} does not have sequence information".format(row["qorf_id"]), file = sys.stderr)
            else:
                # dna sequence extraction
                qseq_dna=extract_sequence(orf2fna, row["qorf_id"], start=row["qostart_dna"], end=row["qoend_dna"])
                sseq_dna=extract_sequence(orf2fna, row["sorf_id"], start=row["sostart_dna"], end=row["soend_dna"])
                
                f.write(">{}:{}\n".format(key, "qseq_dna"))
                f.write(str(qseq_dna)+'\n')
                f.write(">{}:{}\n".format(key, "sseq_dna"))
                f.write(str(sseq_dna)+'\n')         
            
                # protein sequence extraction
                qseq_pro=extract_sequence(orf2faa, row["qorf_id"], start=row["qostart_pro"], end=row["qoend_pro"])
                sseq_pro=extract_sequence(orf2faa, row["sorf_id"], start=row["sostart_pro"], end=row["soend_pro"])
                qseq_pro=str(qseq_pro).replace('*','X')
                sseq_pro=str(sseq_pro).replace('*','X')
                
                # 2. output extracted sequences
                f.write(">{}:{}\n".format(key, "qseq_pro"))
                f.write(str(qseq_pro)+'\n')
                f.write(">{}:{}\n".format(key, "sseq_pro"))
                f.write(str(sseq_pro)+'\n')
    print("DONE: output {}".format(seqFilepath))

def read_fasta(filepath_lst):
    orf2record={}
    for filepath in filepath_lst:
        for record in SeqIO.parse(filepath, "fasta"):
            orfId = record.id
            orf2record[orfId]=record
    return orf2record
    
if __name__=="__main__":
    target=sys.argv[1]
    strain=sys.argv[2]
    annotationType="refseq"

    ovrFilepath="./out/{}/{}_ovr.csv".format(target, strain)
    seqFilepath="./out/{}/{}_ovr.fasta".format(target, strain)
    
    sbjctFilepath="/data/mitsuki/data/denovo/{}/annotation/{}/fna/{}.fna".format(target, annotationType, strain)
    queryFilepath="../blastn/query/{}/{}.fna".format(target, strain)
    orf2fna = read_fasta([sbjctFilepath, queryFilepath])
    orf2faa = read_fasta([sbjctFilepath.replace("fna", "faa"), queryFilepath.replace("fna", "faa")])
    main(target, ovrFilepath, seqFilepath, orf2fna, orf2faa)
