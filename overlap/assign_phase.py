#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62

sys.path.append("../helper")
from phase import PhaseController

def get_six_frame(dnaSeq):
    """
    get 6 frame translation of qseq_dna
    """
    qtrans_lst=[]
    for strand in [1,-1]:
        for gap in range(3):
            start, end = gap, len(dnaSeq)
            while end%3 != gap:
                end -= 1
            if strand == 1:
                subSeq = dnaSeq[gap:end]
            elif strand == -1 :
                subSeq = dnaSeq[gap:end].reverse_complement()
            qtrans_lst.append(subSeq.translate(table=11))
    return qtrans_lst

def get_phase_dct(overlapId, id2rec):
    """
    assign relative reading frame between the 2 sequence of the overlap, by caluculating protein alignment score against 6 frame translation
    """
    
    # extract sequence information
    qseq_dna=id2rec["{}:qseq_dna".format(overlapId)].seq
    sseq_dna=id2rec["{}:sseq_dna".format(overlapId)].seq # never been used for now
    qseq_pro=id2rec["{}:qseq_pro".format(overlapId)].seq
    sseq_pro=id2rec["{}:sseq_pro".format(overlapId)].seq

    # calcurate alignment score against 6 frame translation of qseq_dna (can be sseq_dna though)
    qtrans_lst=get_six_frame(qseq_dna)
    qscore_lst, sscore_lst=[], []
    for qtrans in qtrans_lst:
        qtrans = str(qtrans).replace('*','X')
        alns_pro=pairwise2.align.globalds(qtrans, qseq_pro, matrix, -10, -0.5)
        qscore_lst.append(alns_pro[0][2])
        alns_pro=pairwise2.align.globalds(qtrans, sseq_pro, matrix, -10, -0.5)
        sscore_lst.append(alns_pro[0][2])

    #assign reading frame based on max alignment score
    qrf = np.argmax(qscore_lst)
    srf = np.argmax(sscore_lst)
    
    #update dct and return
    dct={"overlap_id" : overlapId}
    for i in range(6):
        dct["qscore{}".format(i)]=qscore_lst[i]
        dct["sscore{}".format(i)]=sscore_lst[i]
    pc=PhaseController()
    dct["relative"]=pc.phase_lst[pc.relative_int(qrf, srf)]
    return dct

def main(seqFilepath, ovrFilepath, phaseFilepath):
    id2rec={}
    for seqRec in SeqIO.parse(seqFilepath, "fasta"):
        id2rec[seqRec.id]=seqRec
    print("DONE: load {} sequence information".format(int(len(id2rec)/4)))
    
    ovr_df=pd.read_csv(ovrFilepath)
    high_df=ovr_df[(ovr_df["score_pro"] < 20) & (ovr_df["score_dna"] > 100)].copy()  # filter high quality
    print("DONE: filter {}/{} promissing hit.".format(high_df.shape[0], ovr_df.shape[0]))
    
    dct_lst=[]
    for overlapId in high_df["overlap_id"]:
        dct_lst.append(get_phase_dct(overlapId, id2rec))
    phase_df=pd.DataFrame(dct_lst)
    column_lst=["overlap_id", "relative"] + ["qscore{}".format(i) for i in range(6)] + ["sscore{}".format(i) for i in range(6)]
    phase_df = phase_df[column_lst]
    phase_df.to_csv(phaseFilepath, index=False)
    print("DONE: output  {}".format(phaseFilepath))
    
    ovr_df = pd.merge(ovr_df, phase_df[["overlap_id", "relative"]], on = "overlap_id", how = "left")
    ovr_df.to_csv(ovrFilepath , index=False)
    print("DONE: update {}".format(ovrFilepath))

if __name__=="__main__":
    target=sys.argv[1]
    strain=sys.argv[2]
    
    direc="./out/{}".format(target)
    seqFilepath="{}/{}_ovr.fasta".format(direc, strain)
    ovrFilepath="{}/{}_ovr.csv".format(direc, strain)
    phaseFilepath="{}/{}_phase.csv".format(direc, strain)
    main(seqFilepath, ovrFilepath, phaseFilepath)