#!/usr/bin/env python3

import sys
import pandas as pd
from ScoreDbController import ScoreDbController

sys.path.append("../helper")
from gff import read_gff

def main(target, strain_lst):
    print("START: query about {} strains".format(len(strain_lst)))
    dct_lst=[]
    for strain in strain_lst:
        refseqFilepath="/data/mitsuki/data/denovo/{}/annotation/refseq/gff/{}.gff".format(target, strain)
        refseq_df=read_gff(refseqFilepath, ["orf_id"])
        dbFilepath="/data/mitsuki/data/denovo/{}/annotation/prodigal/sup/{}.sq3".format(target, strain)
        sdc=ScoreDbController(dbFilepath)
        
        infoCount = 0
        for _, row in refseq_df.iterrows():
            dct = {"orf_id": row["orf_id"]}
            info_dct = sdc.info(row["seqname"], row["start"], row["end"])
            if len(info_dct) > 0:
                infoCount += 1
                dct.update(info_dct)
            dct_lst.append(dct)
        print("\tDONE: found information for {}/{} genes in {}".format(infoCount, refseq_df.shape[0], strain))
    score_df=pd.DataFrame(dct_lst)
    score_df=score_df[["orf_id", "Beg", "End", "Std", "Total", "CodPot", "StrtSc", "Codon", "RBSMot", "Spacer", "RBSScr", "UpsScr", "TypeScr", "GCCont"]]
    
    outFilepath="../data/{}/orf2score.csv".format(target)
    score_df.to_csv(outFilepath, index=False)
    print("DONE: output {}".format(outFilepath))

if __name__=="__main__":
    target = sys.argv[1]
    strainFilepath="../data/{}/strain.lst".format(target)
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    main(target, strain_lst)
