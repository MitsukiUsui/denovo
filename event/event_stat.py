#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from ete3 import Tree

class EventStat:
    """
    class to calcurate distance between to ortholog family
    EventStat is not appropriate name then
    """
    
    def __init__(self, target):
        self.strain_lst = self.get_strain_lst(target)
        self.cluster_df = self.get_cluster_df(target)
        self.distance_mat = self.get_distance_mat(target)
        
    def get_strain_lst(self, target):
        strainFilepath="../data/{}/strain.lst".format(target)
        strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
        return strain_lst
    
    def get_cluster_df(self, target):
        clusterFilepath="../data/{}/cluster.tsv".format(target)
        cluster_df=pd.read_csv(clusterFilepath, sep="\t", dtype="object")
        return cluster_df
    
    def get_distance_mat(self, target):
        distFilepath = "../data/{}/distance.npy".format(target)
        try:
            distance_mat = np.load(distFilepath)
            print("DONE: load distance matrix from {}".format(distFilepath))
        except FileNotFoundError:
            distance_mat = self.calc_distance_mat(target)
            np.save(distFilepath, distance_mat)
            print("DONE: calc & save distance matrix to {}".format(distFilepath))
            
        for i in range(distance_mat.shape[0]): #to exclude itself in minimum distance calcuration in family_distance
            distance_mat[i][i] = 9999
        return distance_mat
    
    def calc_distance_mat(self, target):
        phbFilepath="../data/{}/cluster.phb".format(target)
        t=Tree(phbFilepath)

        print("START: calc distance matrix between {} strain".format(len(self.strain_lst)))
        distance_mat=-np.ones((len(self.strain_lst), len(self.strain_lst)))
        for i, node1 in enumerate(self.strain_lst):
            for j, node2 in enumerate(self.strain_lst):
                if i != j:
                    distance_mat[i,j]=t.get_distance(node1, node2)              
        return distance_mat
    
    def family_distance(self,  family1, family2):
        try:
            msk1 = self.cluster_df[self.cluster_df["family"] == family1][self.strain_lst].iloc[0].notnull().values
            msk2 = self.cluster_df[self.cluster_df["family"] == family2][self.strain_lst].iloc[0].notnull().values
            return self.distance_mat[msk1, :][:, msk2].min()
        except IndexError:
            print("ERROR: invalid family name {}, {}".format(family1, family2))
            return -1
            

if __name__=="__main__":
    target = sys.argv[1]
    family1 = sys.argv[2]
    family2 = sys.argv[3]
    eventStat = EventStat(target)
    print(eventStat.family_distance(family1, family2))
