#!/usr/bin/env python3

"""
select queryies for mmseqs search against nr
"""

import sys
import os
from Bio import SeqIO
from collections import defaultdict

sys.path.append("../helper")
from myio import get_strain_lst, get_cluster_df

target = sys.argv[1]
strain_lst = get_strain_lst(target)
cluster_df = get_cluster_df(target)

inDirec = "/data/mitsuki/data/denovo/{}/annotation/refseq/faa".format(target)
outDirec = "/data/mitsuki/out/altorf/denovo/trg/{}/mmseqs".format(target)

os.makedirs(outDirec, exist_ok=True)
outFilepath = "{}/query.faa".format(outDirec)

print("START: select representative query randomly for each {} ortholog families".format(cluster_df.shape[0]), flush=True)
query_dct = defaultdict(list) #key: strain, val: list of orfs
for _, row in cluster_df.iterrows():
    sample = row[strain_lst].dropna().sample(n=1)
    strain = sample.index[0]
    orfid = sample.values[0].split(",")[0]
    query_dct[strain].append(orfid)

query_lst = []
for strain in strain_lst:
    orfid2rec = {}
    faaFilepath = "/data/mitsuki/data/denovo/{}/annotation/refseq/faa/{}.faa".format(target, strain)
    for rec in SeqIO.parse(faaFilepath, "fasta"):
        orfid2rec[rec.id] = rec
    for orfid in query_dct[strain]:
        query_lst.append(orfid2rec[orfid])
assert len(query_lst) == cluster_df.shape[0]

with open(outFilepath, 'w') as f:
    for query in query_lst:
        SeqIO.write(query, f, "fasta")
print("DONE: output {} queries to {}".format(len(query_lst), outFilepath), flush=True)
