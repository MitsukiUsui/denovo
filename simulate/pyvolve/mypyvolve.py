#!/usr/bin/env python3

import pyvolve
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import ete3

sys.path.append("../../helper")
from gff import read_gff

def get_evolved(evolver):
    outFilepath="simulated_alignment.fasta"
    evolver(seqfile=outFilepath, ratefile=None, infofile=None)
    record_lst = []
    for record in SeqIO.parse(outFilepath, "fasta"):
        record_lst.append(record)
    return record_lst

def main(strain, seedFilepath, gffFilepath):
    for record in SeqIO.parse(seedFilepath, "fasta"):
        seedRec = record
        break
    gff_df=read_gff(gffFilepath)
    
    #get all the shuffle region
    prv =0
    pos_lst = []
    for _, row in gff_df.iterrows():
        pos_lst.append(("nc", prv, row["start"] - 1, "+"))
        pos_lst.append(("c", row["start"] - 1, row["end"], row["strand"])) 
        prv = row["end"]
    pos_lst.append(("nc", prv, len(seedRec), "+"))
   

    # configuration for evolution
    treeFilepath="tmp.tree"
    mytree = pyvolve.read_tree(file = treeFilepath)
    ncm = pyvolve.Model("nucleotide") # non-coding model
    cm= pyvolve.Model("ECMrest") # coding model

    
    outputSeq_lst=[Seq("") for _ in range(4)] # assuming tree has 4 nodes
    for pos in pos_lst:
        category, start, end, strand  = pos

        # get rootSeq according to start, end, strand info
        rootSeq = seedRec.seq[start:end]
        if strand == "-":
            rootSeq = rootSeq.reverse_complement()
        rootSeq = str(rootSeq)

        # get simulated sequences
        if category == "nc":
#            partition = pyvolve.Partition(models = ncm, root_sequence = rootSeq)
#            evolver = pyvolve.Evolver(partition = partition, tree = mytree)
#            rec_lst = get_evolved(evolver)
            rec_lst = [SeqRecord(Seq(rootSeq)) for _ in range(4)]
        elif category == "c":
            partition = pyvolve.Partition(models = cm, root_sequence = rootSeq[3:-3]) #remove start & stop codon
            evolver = pyvolve.Evolver(partition = partition, tree = mytree)
            rec_lst = get_evolved(evolver)
            for rec in rec_lst:
                rec.seq = rootSeq[:3] + rec.seq + rootSeq[-3:]  #add last stop codon back
        assert len(rec_lst) == len(outputSeq_lst)

        # concat to outputSeq_lst
        for i, rec in enumerate(rec_lst):
            simSeq = rec.seq
            if strand == "-":
                simSeq = simSeq.reverse_complement()
            outputSeq_lst[i] += simSeq
    
    for i, outputSeq in enumerate(outputSeq_lst):
        genomeId = "{}_sim{}".format(strain, i+1)
        outFilepath = "../data/dnaseq/{}.dnaseq".format(genomeId)
        with open(outFilepath, "w") as f:
            seqname = "{}:seq".format(genomeId)
            rec = SeqRecord(outputSeq, id = seqname, description ="")
            SeqIO.write(rec, f, "fasta")
        print("DONE: output {}".format(outFilepath))

if __name__=="__main__":
    strain=sys.argv[1]
    seedFilepath="../data/dnaseq/{}_seed.dnaseq".format(strain)
    gffFilepath=gffFilepath="../data/gff/{}_seed.gff".format(strain)
    main(strain, seedFilepath, gffFilepath)