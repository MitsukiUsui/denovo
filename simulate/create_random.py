#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random

def codon_lookup():
    bases = ["T", "C", "A", "G"]
    id2codon = {}
    codon2id = {}
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codonId = i * 16 + j * 4 + k
                codon = bases[i] + bases[j] + bases[k]
                id2codon[codonId] = codon
                codon2id[codon] = codonId
    return id2codon, codon2id

def revcomp(s):
    return str(Seq(s).reverse_complement())

def get_prv2oks(strand, phase):
    """
    given previous codon-id, return a list of codon-id, which is ok to append after without introducing stop codons in +0 or given strand-phase
    """
    
    assert strand in ("+", "-")
    assert phase in (0, 1, 2)
    
    id2codon, codon2id = codon_lookup()
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    stops = [codon2id[codon] for codon in table.stop_codons]
    if strand == "+":
        pstops = stops # p for phased 
    else:
        pstops = [codon2id[revcomp(codon)] for codon in table.stop_codons]
    
    prv2oks = {}
    if phase == 0:
        prv2oks[-1] = [idx for idx in range(64) if idx not in stops + pstops]
    else:
        prv2oks[-1] = [idx for idx in range(64) if idx not in stops] #-1 for first element
    
    for prv in range(64):
        oks = []
        for idx in range(64):
            codon = id2codon[idx]
            if phase == 0:
                pcodon = codon
            else:
                pcodon = id2codon[prv][(phase - 3):] + codon[:phase]
            pidx = codon2id[pcodon]
            if (idx not in stops) and (pidx not in pstops): # if both codon & phased-codon is not in the stop list respectively
                oks.append(idx)
        prv2oks[prv] = oks
    return prv2oks

def main(aaLength = 100, seqCount = 100):
    print("START: create 5 x {} simulated sequence. ({} aa)".format(seqCount, aaLength))
    id2codon, codon2id = codon_lookup()
    
    for strand, phase in [("+", 0), ("+", 1), ("+", 2), ("-", 0), ("-", 1), ("-", 2)]:
        prv2oks = get_prv2oks(strand, phase)
        seqRec_lst = []
        for simId in range(seqCount): #create total 100 random sequences 
            prv = -1
            seq = ""
            for _ in range(aaLength):
                oks = prv2oks[prv]
                nxt = random.choice(oks)
                seq += id2codon[nxt]
                prv = nxt

            seqname = "sim{}{}_{}".format(strand, phase, simId)
            seqRec = SeqRecord(Seq(seq), id = seqname, description = "")
            seqRec_lst.append(seqRec)
     
        outFilepath = "./out/sim{}{}.fna".format(strand, phase)
        with open(outFilepath, "w") as f:
            for seqRec in seqRec_lst:
                SeqIO.write(seqRec, f, "fasta")
        print("DONE: output {} seqs to {}".format(seqCount, outFilepath))

if __name__=="__main__":
    main()
    
