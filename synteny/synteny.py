#!/usr/bin/env python3

import pandas as pd
import sys
from collections import Counter

sys.path.append("../helper")
from gff import read_gff

def calc_synteny(neighbor_lst):
    c=Counter(neighbor_lst)
    top = c.most_common()[:6]  #select 6 abundant family in the list
    
    neighbor = ",".join([e[0] for e in top])
    synteny = sum([e[1] for e in top]) / len(neighbor_lst)
    return neighbor, synteny

def calc_intra_synteny(family1, family2, syn_df):
    neighbor1_lst = syn_df[syn_df["family"]==family1]["neighbor"].iloc[0].split(",")
    neighbor2_lst = syn_df[syn_df["family"]==family2]["neighbor"].iloc[0].split(",")
    c = Counter(neighbor1_lst+neighbor2_lst)
    commonCount=0
    for k, v in c.items():
        if v==1:
            pass
        elif v==2:
            commonCount+=2
        else:
            print("ERROR: something wrong with {} of {}".format(family1, family2))
            print("\t{}\n\t{}".format(",".join(neighbor1_lst), ",".join(neighbor2_lst)))
            return 0
    return commonCount/len(neighbor1_lst+neighbor2_lst)

def get_orf2synteny(strain_lst, gffDirec):
    orf2synteny = {}
    for strain in strain_lst:
        gffFilepath="{}/{}.gff".format(gffDirec, strain)
        gff_df=read_gff(gffFilepath, ["orf_id", "family"])
        
        for seqname in set(gff_df["seqname"]):
            filtered_df = gff_df[gff_df["seqname"]==seqname].copy()
            filtered_df = filtered_df.sort_values(by=["start", "end"])
            filtered_df = filtered_df.reset_index(drop=True)

            family_lst = list(filtered_df["family"])
            for key, row in filtered_df.iterrows():
                pre_lst = family_lst[max(0, key - 3): key]
                post_lst = family_lst[key + 1: key + 4] 

                dct = {}
                dct["neighbor"] = pre_lst + post_lst
                dct["left"], dct["right"] = "", ""

                if len(pre_lst) > 0:
                    if row["strand"] == "+":
                        dct["left"] = pre_lst[-1]
                    elif row["strand"] == "-":
                        dct["right"] = pre_lst[-1]
                    else:
                        print("ERROR: unknown strand {}".format(row["strand"]))
                if len(post_lst) > 0:
                    if row["strand"] == "+":
                        dct["right"] = post_lst[0]
                    elif row["strand"] == "-":
                        dct["left"] = post_lst[0]
                    else:
                        print("ERROR: unknown strand {}".format(row["strand"]))
                orf2synteny[row["orf_id"]] = dct    
    return orf2synteny 
    
def get_synteny_df(cluster_df, orf2synteny, strain_lst):
    dct_lst=[]
    for _, row in cluster_df.iterrows():
        left_lst, right_lst, neighbor_lst = [], [], []
        for orfIds in row[strain_lst].dropna():
            for orfId in orfIds.split(","):
                left_lst.append(orf2synteny[orfId]["left"])
                right_lst.append(orf2synteny[orfId]["right"])
                neighbor_lst += orf2synteny[orfId]["neighbor"]

        dct = {"family": row["family"]}
        dct["left"], dct["left_per"] = Counter(left_lst).most_common()[0]
        dct["left_per"] /= len(left_lst)
        dct["right"], dct["right_per"] = Counter(right_lst).most_common()[0]
        dct["right_per"] /= len(right_lst)
        
        try:
            dct["neighbor"], dct["synteny"]=calc_synteny(neighbor_lst)
        except:
            print("ERROR: {}".format(row["family"]))
        dct_lst.append(dct)
    syn_df = pd.DataFrame(dct_lst)
    syn_df=syn_df[["family", "left_per", "right_per", "synteny", "left", "right", "neighbor"]]
    return syn_df
    

def main(strainFilepath, clusterFilepath, outFilepath):
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    cluster_df=pd.read_csv(clusterFilepath, delimiter="\t")
    
    gffDirec="/data/mitsuki/data/ortho/bacillus/gff"
    print("START: gether synteny information for {} strains".format(len(strain_lst)))
    orf2synteny=get_orf2synteny(strain_lst, gffDirec)
    print("START: calculate synteny score for {} families".format(cluster_df.shape[0]))
    syn_df = get_synteny_df(cluster_df, orf2synteny, strain_lst)
    syn_df.to_csv(outFilepath, sep="\t", index=False)
    print("DONE: output to {}".format(outFilepath))

if __name__=="__main__":
    target="bacillus"
    strainFilepath="../data/{}/strain.lst".format(target)
    clusterFilepath="../data/{}/cluster.tsv".format(target)
    outFilepath = "../data/{}/synteny.tsv".format(target)
    main(strainFilepath, clusterFilepath, outFilepath)
