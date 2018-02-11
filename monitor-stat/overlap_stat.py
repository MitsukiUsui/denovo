#!/usr/bin/env python3

import sys
import pandas as pd
from collections import defaultdict

sys.path.append("../helper")
from myio import get_strain_lst
from interval import *

phase_lst = ["+0", "+1", "+2", "-0", "-1", "-2"]

def main(target, outFilepath):
    strain_lst = get_strain_lst(target)
    print("START: collect overlap stat for {} strains".format(len(strain_lst)), flush = True)
    
    dct_lst = []
    for strain in strain_lst:
        dct = {"strain": strain}
        ovrFilepath = "../overlap/out/{}/{}_ovr.csv".format(target, strain)
        ovr_df = pd.read_csv(ovrFilepath, dtype = {"relative":object})
        
        count_dct = {}
        for phase in phase_lst + ["undefined", "all"]:
            count_dct[phase] = 0
            
        chrName_lst = list(set(ovr_df["chr_name"]))
        for chrName in chrName_lst:
            filter_df = ovr_df[(ovr_df["sfamily"].notnull()) & (ovr_df["chr_name"] == chrName)] #filter pseudogene 

            phase2interval = defaultdict(list) # key: phase, value: list of interval
            for _, row in filter_df.iterrows():
                if isinstance(row["relative"], str):
                    phase2interval[row["relative"]].append(Interval(row["ostart"], row["oend"]))
                else:
                    assert row["olength"] <= 10
                    phase2interval["undefined"].append(Interval(row["ostart"], row["oend"]))

            # update count_dct
            all_lst = []
            for phase in phase_lst + ["undefined"]:
                all_lst += phase2interval[phase]
                count_dct[phase] += interval_sum(phase2interval[phase])
            count_dct["all"] += interval_sum(all_lst)
        
        dct.update(count_dct)
        dct_lst.append(dct)

    stat_df = pd.DataFrame(dct_lst)
    stat_df = stat_df[["strain", "all", "+0", "+1", "+2", "-0", "-1", "-2", "undefined"]]
    stat_df.to_csv(outFilepath, index = False)
    print("DONE: output {}".format(outFilepath), flush = True)
   
    
if __name__ == "__main__":
    target = sys.argv[1]
    outFilepath = "./out/overlap/{}.csv".format(target)
    main(target, outFilepath)