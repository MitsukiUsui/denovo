#!/usr/bin/env python3

import sys
import pandas as pd

sys.path.append("../helper")
from myio import *

target = sys.argv[1]
strain_lst = get_strain_lst(target)

# check result
for strain in strain_lst:
    hitFilepath = "./result/{}/{}.csv".format(target, strain)
    try:
        hit_df=pd.read_csv(hitFilepath)
    except FileNotFoundError:
        print("ERROR: {} does not exist".format(hitFilepath), file = sys.stderr)
        exit(1)
exit(0)
