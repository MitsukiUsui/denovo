#!/usr/bin/env python3

import sys
import pandas as pd

sys.path.append("../helper")
from myio import *

target=sys.argv[1]
strain_lst = get_strain_lst(target)

# check sq3?

#outFilepath="../data/{}/family2score.csv".format(target)
#try:
#    out_df = pd.read_csv(outFilepath)
#except FileNotFoundError:
#    print("ERROR: {} dose not exist".format(eventFilepath), file = sys.stderr)
#    exit(1)
exit(0)
    
