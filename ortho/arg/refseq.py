#!/usr/bin/env python3

import pandas as pd
import sys
sys.path.append("../helper")
from myutil.myutil import myrun

target=sys.argv[1]
catalogFilepath="../data/{}/catalog.tsv".format(target)
catalog_df=pd.read_csv(catalogFilepath, sep="\t")
for _, row in catalog_df.iterrows():
    basename = row["ftp_path"].split("/")[-1]
    genomeId = row["genome_id"]
    cmd = "{},{},{}".format(target, basename, genomeId)
    print(cmd)
#    myrun(cmd)
