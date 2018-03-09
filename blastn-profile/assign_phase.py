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
    qseq_pro = str(qseq_pro).replace('U', 'X')
    sseq_pro = str(sseq_pro).replace('U', 'X')

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

def calc_score_dna(qseq_dna, sseq_dna, strand):
    if strand == 1:
        alns_dna=pairwise2.align.globalms(qseq_dna, sseq_dna, 2, -1, -.5, -.1)
    elif strand == -1:
        alns_dna=pairwise2.align.globalms(qseq_dna, sseq_dna.reverse_complement(), 2, -1, -.5, -.1)
    return alns_dna[0][2]

def main(ovrFilepath, seqFilepath, phaseFilepath):
    
    ovr_df = pd.read_csv(ovrFilepath)
    id2seq = {}
    for rec in SeqIO.parse(seqFilepath, "fasta"):
        id2seq[rec.id] = rec.seq
    
    print("START: assign phase to {} overlaps".format(ovr_df.shape[0]))
    dct_lst = []
    for _, row in ovr_df.iterrows():
        if row["olength"] > 10:
            qseq_dna = id2seq["{}:qseq_dna".format(row["region_id"])]
            sseq_dna = id2seq["{}:sseq_dna".format(row["region_id"])]
            qseq_pro = id2seq["{}:qseq_pro".format(row["region_id"])]
            sseq_pro = id2seq["{}:sseq_pro".format(row["region_id"])]

            dct = {"region_id": row["region_id"],
                    "score_dna": calc_score_dna(qseq_dna, sseq_dna, row["qstrand"] * row["sstrand"])}
            dct.update(get_phase_dct(qseq_dna, sseq_dna, qseq_pro, sseq_pro))
            dct_lst.append(dct)

    phase_df=pd.DataFrame(dct_lst)
    column_lst=["region_id", "score_dna", "relative"] + ["qscore{}".format(i) for i in range(6)] + ["sscore{}".format(i) for i in range(6)]
    phase_df = phase_df[column_lst]
    phase_df.to_csv(phaseFilepath, index=False)
    print("DONE: output  {}".format(phaseFilepath))

if __name__=="__main__":
    target=sys.argv[1]
    strain=sys.argv[2]
    
    direc = "./out/{}/{}".format(target, strain)
    ovrFilepath = "{}/extract.csv".format(direc)
    seqFilepath = "{}/extract.fasta".format(direc)
    phaseFilepath = "{}/phase.csv".format(direc)
    main(ovrFilepath, seqFilepath, phaseFilepath)
