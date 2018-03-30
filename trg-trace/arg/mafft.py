#!/usr/bin/env python3

import sys
import pandas as pd

sys.path.append("../helper")
from myio import get_strain_lst

target = sys.argv[1]
strain_lst = get_strain_lst(target)
for strain in strain_lst:
    traceFilepath = "./out/{}/{}.trace".format(target, strain)
    trace_df = pd.read_csv(traceFilepath)
    for query in trace_df[trace_df["traceable"]==1]["query"]:
        fnaFilepath = "./align/{}/{}/{}.fna".format(target, strain, query)
        alignFilepath = fnaFilepath.replace(".fna", ".align")
        print("{},{}".format(fnaFilepath, alignFilepath))

