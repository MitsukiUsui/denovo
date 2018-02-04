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

def get_phase_dct(qseq_dna, sseq_dna, qseq_pro, sseq_pro):
    """
    assign relative reading frame between the 2 sequence of the overlap, by caluculating protein alignment score against 6 frame translation of qseq_dna
    """
    
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
    dct = {}
    for i in range(6):
        dct["qscore{}".format(i)]=qscore_lst[i]
        dct["sscore{}".format(i)]=sscore_lst[i]
    pc=PhaseController()
    dct["relative"]=pc.phase_lst[pc.relative_int(qrf, srf)]
    return dct

def main(ovrFilepath, seqFilepath, phaseFilepath):
    rec_lst = []
    for rec in SeqIO.parse(seqFilepath, "fasta"):
        rec_lst.append(rec)
    
    assert len(rec_lst) % 4 == 0
    numOvr = int(len(rec_lst)/4)
    print("START: assign phase to {} overlaps".format(numOvr))
    dct_lst = []
    for i in range(numOvr):
        qseq_dna = rec_lst[4 * i + 0].seq
        sseq_dna = rec_lst[4 * i + 1].seq
        qseq_pro = rec_lst[4 * i + 2].seq
        sseq_pro = rec_lst[4 * i + 3].seq
       
        overlapId = int(rec_lst[4 * i + 0].id.split(':')[0])
        dct = {"overlap_id": overlapId}
        dct.update(get_phase_dct(qseq_dna, sseq_dna, qseq_pro, sseq_pro))
        dct_lst.append(dct)
    phase_df=pd.DataFrame(dct_lst)
    column_lst=["overlap_id", "relative"] + ["qscore{}".format(i) for i in range(6)] + ["sscore{}".format(i) for i in range(6)]
    phase_df = phase_df[column_lst]
    phase_df.to_csv(phaseFilepath, index=False)
    print("DONE: output  {}".format(phaseFilepath))
    
    
    ovr_df = pd.read_csv(ovrFilepath)
    ovr_df = pd.merge(ovr_df, phase_df[["overlap_id", "relative"]], on = "overlap_id", how = "left")
    ovr_df.to_csv(ovrFilepath, index = False)
    print("DONE: append relative column to {}".format(ovrFilepath))

if __name__=="__main__":
    target=sys.argv[1]
    strain=sys.argv[2]
    
    direc="./out/{}".format(target)
    ovrFilepath="{}/{}_ovr.csv".format(direc, strain)
    seqFilepath="{}/{}_ovr.fasta".format(direc, strain)
    phaseFilepath="{}/{}_phase.csv".format(direc, strain)
    main(ovrFilepath, seqFilepath, phaseFilepath)