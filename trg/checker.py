#!/usr/bin/env python3

import sys
import pandas as pd

target=sys.argv[1]
trgFilepath = "/data/mitsuki/out/altorf/denovo/trg/{}/trg.csv".format(target)
try:
    trg_df = pd.read_csv(trgFilepath, dtype={"genus": int})
except FileNotFoundError:
    print("ERROR: {} dose not exist".format(trgFilepath), file = sys.stderr)
    exit(1)
except KeyError:
    print("ERROR: {} dose not column {}".format(trgFilepath), file = sys.stderr)
    exit(1)
exit(0)
    
