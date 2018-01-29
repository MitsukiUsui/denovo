#!/usr/bin/env python3

import sys
import pandas as pd
from collections import defaultdict, Counter
import numpy as np

import sys
sys.path.append("../synteny")
from synteny import calc_inter_synteny
sys.path.append("../helper")
from phase import PhaseController

def get_strain_lst(target):
    strainFilepath="../data/{}/strain.lst".format(target)
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    return strain_lst

def get_cluster_df(target):
    clusterFilepath="../data/{}/cluster.tsv".format(target)
    cluster_df=pd.read_csv(clusterFilepath, sep="\t", dtype="object")
    return cluster_df

def get_synteny_df(target):
    syntenyFilepath="../data/{}/synteny.tsv".format(target)
    synteny_df=pd.read_csv(syntenyFilepath, sep = "\t")
    return synteny_df

def get_all_df(target, strain_lst):
    dct_lst=[]
    for strain in strain_lst:
        filepath="./out/{}/{}_ovr.csv".format(target, strain)
        try:
            ovr_df=pd.read_csv(filepath, dtype = {"relative":object})
            ovr_df["sstrain"]=strain
            dct_lst+=ovr_df.to_dict("records")
        except FileNotFoundError:
            print("WARN: {} does not exist".format(filepath))
    all_df=pd.DataFrame(dct_lst)
    column_lst=["overlap_id", "sstrain", "qstrain",  "olength", "score_dna", "score_pro", 
                     "qfamily", "sfamily", "qorf_id", "sorf_id", "relative"]
    all_df=all_df[column_lst]
    return all_df

def get_family2lineage(cluster_df):
    family2lineage={}
    for _, row in cluster_df.iterrows():
        family2lineage[row["family"]]=int(row["lineage"])
    return family2lineage

class Stat:
    def __init__(self):
        self.rows=[]
        self.df = None

    def update(self, row):
        if len(self.rows) == 0:
            self.left = min(row["qfamily"], row["sfamily"])
            self.right = max(row["qfamily"], row["sfamily"])
        else:
            assert row["qfamily"] in (self.left, self.right)
            assert row["sfamily"] in (self.left, self.right)
        self.rows.append(row)
        
    def calc(self):
        if self.df is None:
            self.df = pd.DataFrame(self.rows)
            if self.df.shape[0] > 0:
                self.ldf = self.df[self.df["qfamily"] == self.left]
                self.rdf = self.df[self.df["qfamily"] == self.right]
            
    def get_hitCount(self):
        self.calc()
        if self.df.shape[0] == 0:
            return 0, 0
        else:
            ret1 = self.ldf.shape[0]
            ret2 = self.rdf.shape[0]
            return ret1, ret2
        
    def get_sbjctCount(self):
        self.calc()
        if self.df.shape[0] == 0:
            return 0, 0
        else:
            ret1 = len(set(self.ldf["sstrain"]))
            ret2 = len(set(self.rdf["sstrain"]))
            return ret1, ret2
        
    def get_aveLength(self):
        self.calc()
        if self.df.shape[0] == 0:
            return 0, 0
        else:
            ret1 = self.ldf["olength"].mean()
            ret2 = self.rdf["olength"].mean()
            return ret1, ret2
        
    def get_phase(self):
        self.calc()
        if self.df.shape[0] == 0:
            return None, None 
        else:
            pc = PhaseController()
            phase_lst=[]
            for phase in self.ldf["relative"]:
                phase_lst.append(phase)
            for phase in self.rdf["relative"]:
                phase_lst.append(pc.revops(phase))

            phase, confidence  = Counter(phase_lst).most_common()[0]
            confidence /= len(phase_lst)
            return phase, confidence

def main(target, outFilepath):
    strain_lst=get_strain_lst(target)
    all_df=get_all_df(target, strain_lst)
    high_df=all_df[(all_df["score_pro"] < 20) & (all_df["score_dna"] > 100)].copy()  # filter high quality
    print("DONE: filter {}/{} promissing hit.".format(high_df.shape[0], all_df.shape[0]))
    cluster_df=get_cluster_df(target)
    synteny_df=get_synteny_df(target)
    family2lineage=get_family2lineage(cluster_df)
    
    # update stat object
    key2stat = defaultdict(Stat)
    for _, row in high_df.iterrows():
        key = (min(row["qfamily"], row["sfamily"]), max(row["qfamily"], row["sfamily"]))
        key2stat[key].update(row)
    
    dct_lst=[]
    for key, stat in key2stat.items():
        assert key[0] < key[1]
        dct = {"left": key[0], "right": key[1]}
        dct["left_hit_count"], dct["right_hit_count"] = stat.get_hitCount()
        dct["left_sbjct_count"], dct["right_sbjct_count"] = stat.get_sbjctCount()
        dct["left_average_length"], dct["right_average_length"] = stat.get_aveLength()
        dct["phase"], dct["phase_conf"] = stat.get_phase()
        try:
            dct["inter_synteny"] = calc_inter_synteny(key[0], key[1], synteny_df)
        except:
            dct["inter_synteny"] = -1
        dct["left_lineage"] = family2lineage[key[0]]
        dct["right_lineage"]=family2lineage[key[1]]
        dct["both_lineage"]=(cluster_df[cluster_df["family"].isin(key)][strain_lst].isnull().sum(axis=0) ==0).sum()
        dct_lst.append(dct)
    out_df=pd.DataFrame(dct_lst)
    out_df=out_df[["left", "right", "left_lineage", "right_lineage", "both_lineage", "inter_synteny",
                         "left_hit_count", "right_hit_count", "left_sbjct_count", "right_sbjct_count", 
                         "left_average_length", "right_average_length", "phase", "phase_conf"]]
    out_df.index.name="event_id"
    out_df.to_csv(outFilepath)
    print("DONE: summarize promissing hit into {}".format(outFilepath))

if __name__=="__main__":
    target=sys.argv[1]
    outFilepath="./out/{}/event.csv".format(target)
    main(target, outFilepath)
