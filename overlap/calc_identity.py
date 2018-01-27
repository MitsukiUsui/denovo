#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import subprocess
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo as matlist

def extract_sequence(orf2record, orfId, start, end):
    return orf2record[orfId].seq[start:end]

def main(target, overlapFilepath, logFilepath, orf2fna, orf2faa):
    overlap_df=pd.read_csv(overlapFilepath, index_col=0)
    
    # calc alignment score for geneseq & proteinseq
    matrix = matlist.blosum62
    scoreDna_lst, scorePro_lst = [], []
    print("TOTAL {} alignments".format(overlap_df.shape[0]), flush = True)
    
    with open(logFilepath, 'w') as f:
        for key,row in overlap_df.iterrows():
            if key%100==0:
                print("\t{}".format(key), flush = True)

            if (row["olength"] < 10):
                print("DEBUG: overlap length too short.")
                scoreDna_lst.append(0)
                scorePro_lst.append(0)
            elif (row["qorf_id"] not in orf2fna) or (row["qorf_id"] not in orf2faa) or \
                 (row["sorf_id"] not in orf2fna) or (row["sorf_id"] not in orf2faa):
                print("DEBUG: {} or {} does not have sequence information".format(row["qorf_id"], row["sorf_id"]))
                scoreDna_lst.append(0)
                scorePro_lst.append(0)
            else:
                # dnaseq alignment
                # 1. sequence extraction
                qseq_dna=extract_sequence(orf2fna, row["qorf_id"], start=row["qostart_dna"], end=row["qoend_dna"])
                sseq_dna=extract_sequence(orf2fna, row["sorf_id"], start=row["sostart_dna"], end=row["soend_dna"])

                # 2. output extracted sequences
                f.write(">{}:{}\n".format(key, "qseq_dna"))
                f.write(str(qseq_dna)+'\n')
                f.write(">{}:{}\n".format(key, "sseq_dna"))
                f.write(str(sseq_dna)+'\n')
                        
                # 3. alignment calculation
                try:
                    if row["qstrand"] * row["sstrand"]==1:
                        alns_dna=pairwise2.align.globalms(qseq_dna, sseq_dna, 2, -1, -.5, -.1)
                    elif row["qstrand"] * row["sstrand"]==-1:
                        alns_dna=pairwise2.align.globalms(qseq_dna, sseq_dna.reverse_complement(), 2, -1, -.5, -.1)
                    else:
                        print("ERROR: bad strand format ({}, {}).".format(row["qstrand", "sstrand"]))
                        exit(1)
                    scoreDna_lst.append(alns_dna[0][2])
                except (IndexError, SystemError):
                    print("WARN: failed calc DNA alignment score for overlap_id = {}".format(row["overlap_id"]))
                    scoreDna_lst.append(0)                
            
                #proteinseq alignment
                # 1. sequence extraction
                qseq_pro=extract_sequence(orf2faa, row["qorf_id"], start=row["qostart_pro"], end=row["qoend_pro"])
                sseq_pro=extract_sequence(orf2faa, row["sorf_id"], start=row["sostart_pro"], end=row["soend_pro"])
                qseq_pro=str(qseq_pro).replace('*','X')
                sseq_pro=str(sseq_pro).replace('*','X')
                
                # 2. output extracted sequences
                f.write(">{}:{}\n".format(key, "qseq_pro"))
                f.write(str(qseq_pro)+'\n')
                f.write(">{}:{}\n".format(key, "sseq_pro"))
                f.write(str(sseq_pro)+'\n')
                
                # 3. alignment calculation
                try: 
                    alns_pro=pairwise2.align.globalds(qseq_pro, sseq_pro, matrix, -10, -0.5)
                    scorePro_lst.append(alns_pro[0][2])
                except (IndexError, SystemError):
                    print("WARN: failed calc Protein alignment score for overlap_id = {}".format(row["overlap_id"]))
                    scorePro_lst.append(0)

    overlap_df["score_dna"]=scoreDna_lst
    overlap_df["score_pro"]=scorePro_lst
    overlap_df.to_csv(overlapFilepath)

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
    annotationType="prodigal"

    overlapFilepath="./out/{}/{}_ovr.csv".format(target, strain)
    logFilepath="./out/{}/{}_ovr.fasta".format(target, strain)
    
    sbjctFilepath="/data/mitsuki/data/denovo/{}/annotation/{}/fna/{}.fna".format(target, annotationType, strain)
    queryFilepath="../blastn/query/{}/{}.fna".format(target, strain)
    orf2fna = read_fasta([sbjctFilepath, queryFilepath])
    orf2faa = read_fasta([sbjctFilepath.replace("fna", "faa"), queryFilepath.replace("fna", "faa")])
    main(target, overlapFilepath, logFilepath, orf2fna, orf2faa)
