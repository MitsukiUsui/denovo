#!/usr/bin/env python3

from ete3 import Tree
import pandas as pd
import numpy as np
import re
import sys

def main(clusterFilepath, strainFilepath, phbFilepath, outFilepath):
    cluster_df=pd.read_csv(clusterFilepath, delimiter='\t', dtype="object")
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    t=Tree(phbFilepath)
  
    print("START: calc distance matrix between {} strain".format(len(strain_lst)))
    distance_mat=-np.ones((len(strain_lst), len(strain_lst)))
    for i, node1 in enumerate(strain_lst):
        for j, node2 in enumerate(strain_lst):
            if i!=j:
                distance_mat[i,j]=t.get_distance(node1, node2)

    print("START: create query lookup table".format(len(strain_lst)))
    dct_lst=[]
    for _, row in cluster_df.iterrows():
        lineage=int(row["lineage"])
        if lineage==1 or lineage == len(strain_lst):
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
    direc="../data/{}".format(target)
    
    clusterFilepath="{}/cluster.tsv".format(direc)
    strainFilepath="{}/strain.lst".format(direc)
    phbFilepath="{}/cluster.phb".format(direc)
    outFilepath="{}/query_lookup.csv".format(direc)
    main(clusterFilepath, strainFilepath, phbFilepath, outFilepath)
