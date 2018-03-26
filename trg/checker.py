#!/usr/bin/env python3

import sys
import pandas as pd

target=sys.argv[1]
lcaFilepath = "/data/mitsuki/out/altorf/denovo/trg/{}/lca.csv".format(target)
try:
    lca_df = pd.read_csv(lcaFilepath, dtype={"trg": int})
except FileNotFoundError:
    print("ERROR: {} dose not exist".format(trgFilepath), file = sys.stderr)
    exit(1)
except KeyError:
    print("ERROR: {} dose not column {}".format(trgFilepath), file = sys.stderr)
    exit(1)
exit(0)

