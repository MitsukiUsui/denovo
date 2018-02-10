#!/usr/bin/env python3

from ete3 import Tree
import pandas as pd
import numpy as np
import re
import sys

sys.path.append("../helper")
from myio import *

def main(target, outFilepath):
    strain_lst = get_strain_lst(target)
    cluster_df = get_cluster_df(target)
    distance_mat = get_distance_mat(target)

    print("START: create query lookup table".format(len(strain_lst)))
    dct_lst=[]
    for _, row in cluster_df.iterrows():
        lineage=int(row["lineage"])
        #if lineage==1 or lineage == len(strain_lst):
        if lineage == len(strain_lst):
            pass
        else:
            dct={ "family" : row["family"] }
            msk=row[strain_lst].isnull()
            for sidx in range(len(strain_lst)): #subject index
                if msk[sidx]:
                    x = np.ma.array(distance_mat[sidx], mask=msk)
                    qidx=x.argmin() #query index
                    assert distance_mat[sidx,qidx]>=0
                    orfId=row[strain_lst[qidx]].split(',')[0] #use only first gene as a query 
                    dct[strain_lst[sidx]]=orfId
            dct_lst.append(dct)
    out_df=pd.DataFrame(dct_lst)
    out_df=out_df[["family"]+strain_lst]
    out_df.to_csv(outFilepath, index=False)
    print("DONE: output table to {}".format(outFilepath))

if __name__=="__main__":
    target=sys.argv[1]
    outFilepath="../data/{}/query_lookup.csv".format(target)
    main(target, outFilepath)
