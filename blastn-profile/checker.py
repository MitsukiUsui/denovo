#!/usr/bin/env python3

import sys
import pandas as pd

sys.path.append("../helper")
from myio import get_strain_lst

target=sys.argv[1]
strain_lst = get_strain_lst(target)
for strain in strain_lst:
    fp = "./out/{}/{}/profile.csv".format(target, strain)
    try:
        df = pd.read_csv(fp)
    except FileNotFoundError:
        print("ERROR: {} dose not exist".format(fp), file = sys.stderr)
        exit(1)
exit(0)
    
