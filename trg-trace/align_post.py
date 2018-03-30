#!/usr/bin/env python3

import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

sys.path.append("../helper")
from myio import get_strain_lst

class GapTracker:
    def __init__(self):
        self.gap = {}
        self.gaps = defaultdict(list)

    def start(self, i, name):
        if name not in self.gap.keys():
            self.gap[name] = {"start": i}

    def end(self, i, name):
        if name in self.gap.keys():
            self.gap[name]["end"] = i
            self.gaps[name].append(self.gap[name])
            del self.gap[name]

def parse_align(fp):
    rec_lst = [rec for rec in SeqIO.parse(fp, "fasta")]
    assert len(rec_lst) == 2
    query = rec_lst[0].seq
    sbjct = rec_lst[1].seq
    assert len(query) == len(sbjct)
    return query, sbjct

def is_typical(seq):
    if len(seq)%3 != 0:
        print("WARN: the length is not multiple of 3.", file=sys.stderr)
        return False
    if seq[-3:].translate()!="*":
        print("WARN: the last codon is not stop codon.", file=sys.stderr)
        return False
    if '*' in seq[:-3].translate():
        print("WARN: in-frame stop codon found.", file=sys.stderr)
        return False
    return True

def compare(qseq, sseq):
    """
    compare qseq and sseq
    :ret inframe, indel
    """

    assert len(qseq) == len(sseq)
    length = len(qseq)

    qi = 0 #record query index
    gt = GapTracker()
    inframe = 0

    codon  = {"query": "", "sbjct": ""}
    isValid = True

    # process
    for i in range(length):
        if qseq[i] != '-' and sseq[i] != '-': # if alignment found
            qi += 1
            codon["query"] += qseq[i]
            codon["sbjct"] += sseq[i]
            gt.end(i, "query")
            gt.end(i, "sbjct")
        else:
            isValid = False
            if qseq[i] == '-':
                gt.start(i, "query")
                gt.end(i, "sbjct")
            elif sseq[i] == '-':
                gt.start(i, "sbjct")
                gt.end(i, "query")

        if qi > 0 and qi%3 == 0:
            if isValid:
                assert len(codon["query"]) == 3
                assert len(codon["sbjct"]) == 3
                aa = {"query":  Seq(codon["query"]).translate(),
                      "sbjct": Seq(codon["sbjct"]).translate()}
                if aa["query"]!='*' and aa["sbjct"]=='*':
                    inframe += 1
            codon  = {"query": "", "sbjct": ""}
            isValid = True
    gt.end(length, "query")
    gt.end(length, "sbjct")

    # format indel into 1 string
    indel_lst = []
    for name, prefix in [("query", "-"), ("sbjct", "+")]:
        for gap in gt.gaps[name]:
            if gap["start"]!=0 and  gap["end"]!=length:
                indel_lst.append("{0}{1}({2})".format(prefix, gap["end"]-gap["start"], gap["start"]))
    return inframe, ';'.join(indel_lst)

def main(traceFilepath, alignDirec, outFilepath):
    trace_df = pd.read_csv(traceFilepath)
    query_lst = list(trace_df[trace_df["traceable"]==1]["query"])

    dct_lst = []
    for query in query_lst:
        alignFilepath = "{}/{}.align".format(alignDirec, query)
        qseq, sseq = parse_align(alignFilepath)

        dct = {"query": query}
        if is_typical(Seq(str(qseq).replace('-', ''))):
            inframe, indel = compare(qseq, sseq)
            dct["inframe"] = inframe
            dct["indel"] = indel
        else:
            print("SKIP: atypical CDS {}".format(query))
        dct_lst.append(dct)

    try:
        out_df = pd.DataFrame(dct_lst)
        out_df = out_df[["query", "inframe", "indel"]]
        out_df.to_csv(outFilepath, index=False)
        print("DONE: output {}".format(outFilepath))
    except KeyError:
        pass
if __name__=="__main__":
    target = sys.argv[1]
    strain_lst = get_strain_lst(target)
    for strain in strain_lst:
        traceFilepath = "./out/{}/{}.trace".format(target, strain)
        alignDirec = "./align/{}/{}".format(target, strain)
        outFilepath = "{}/summary.csv".format(alignDirec)
        main(traceFilepath, alignDirec, outFilepath)
