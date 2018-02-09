#!/usr/bin/env python3

import sys
import pandas as pd

target=sys.argv[1]
strainFilepath="../data/{}/strain.lst".format(target)
strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]

# check score
outFilepath="../data/{}/family2score.csv".format(target)
try:
    out_df = pd.read_csv(outFilepath)
except FileNotFoundError:
    print("ERROR: {} dose not exist".format(eventFilepath), file = sys.stderr)
    exit(1)
exit(0)
    
