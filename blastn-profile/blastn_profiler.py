#!/usr/bin/env python3

import sys
import pandas as pd

def main(regionFilepath, phaseFilepath, outFilepath):
    region_df = pd.read_csv(regionFilepath)
    phase_df = pd.read_csv(phaseFilepath, dtype={"relative":object})
    
    join_df = pd.merge(region_df, phase_df[["region_id", "relative"]], on="region_id", how="left")
    assert join_df.shape[0] == region_df.shape[0]
    
    # renew category list
    category_lst = []
    for _, row in join_df.iterrows():
        if row["category"] == "genic":
            if isinstance(row["relative"], str):
                if row["relative"] == "+0":
                    category = "same"
                else:
                    category = "shifted"
            else:
                category = "short"
        else:
            category = row["category"]
        category_lst.append(category)
    join_df["category"] = category_lst
    
    out_df = join_df[["region_id", "hit_id", "qorf_id", "seqname","category", "rqstart", "rqend"]]
    out_df.to_csv(outFilepath, index=False)
    print("DONE: output {}".format(outFilepath))
    
if __name__=="__main__":
    target=sys.argv[1]
    strain=sys.argv[2]

    direc = "./out/{}/{}".format(target, strain)
    #direc = "./out"
    regionFilepath = "{}/region.csv".format(direc)
    phaseFilepath = "{}/phase.csv".format(direc)
    outFilepath = "{}/profile.csv".format(direc)
    main(regionFilepath, phaseFilepath, outFilepath)
