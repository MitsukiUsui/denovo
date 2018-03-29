#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

if __name__=="__main__":
    config_df = pd.read_csv("extract.config")
    outFilepath = "extract.out"

    print("START: output to {}".format(outFilepath))
    with open(outFilepath, "w") as f:
        for _, row in config_df.iterrows():
            seqFilepath = "/data/mitsuki/data/denovo/{}/dnaseq/{}.dnaseq".format(row["target"], row["strain"])

            seq = None
            for rec in SeqIO.parse(seqFilepath, "fasta"):
                if rec.id == row["seqid"]:
                    seq = rec.seq
                    break
            if seq is None:
                print("ERROR: {} not found".format(row["seqid"]))
                exit()

            name = "{}:[{}, {}),{}".format(row["seqid"], row["start"], row["end"], row["reverse"])
            if row["reverse"] == 0:
                outRec = SeqRecord(seq[row["start"]:row["end"]], id = name, description = "")
            else:
                outRec = SeqRecord(seq[row["start"]:row["end"]].reverse_complement(), id = name, description = "")

            SeqIO.write(outRec, f, "fasta")
            print("\t{}".format(name))

