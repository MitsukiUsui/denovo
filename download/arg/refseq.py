#!/usr/bin/env python3

import pandas as pd
import sys
sys.path.append("../helper")

target=sys.argv[1]
catalogFilepath="../data/{}/catalog.tsv".format(target)
catalog_df=pd.read_csv(catalogFilepath, sep="\t")
for _, row in catalog_df.iterrows():
    basename = row["ftp_path"].split("/")[-1]
    genomeid = row["genome_id"]
    print("{},{},{}".format(target, basename, genomeid))
