#!/usr/bin/env python3

import sys
import os
import pandas as pd

sys.path.append("../helper")
from gff import read_gff
from myio import *

target=sys.argv[1]
annotType="refseq"
strain_lst = get_strain_lst(target)

# check cluster.tsv
cluster_df = get_cluster_df(target)
for col in ["family", "lineage", "size"] + strain_lst:
    if not(col in cluster_df.columns):
        print("ERROR: {} dose not have column {}".format(clusterFilepath, col), file = sys.stderr)
        exit(1)

# check gff
for strain in strain_lst:
    gffFilepath = "/data/mitsuki/data/denovo/{}/annotation/{}/gff/{}.gff".format(target, annotType, strain)
    try:
        read_gff(gffFilepath, ["family"])
    except KeyError:
        print("ERROR: {} does not have family information".format(gffFilepath), file = sys.stderr)
        exit(2)
        
exit(0)
