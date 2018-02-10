#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from ete3 import Tree

sys.path.append("../helper")
from myio import *

class EventStat:
    """
    class to calcurate distance between to ortholog family
    EventStat is not appropriate name then
    """
    
    def __init__(self, target):
        self.strain_lst = get_strain_lst(target)
        self.cluster_df = get_cluster_df(target)
        self.distance_mat = get_distance_mat(target)
        for i in range(self.distance_mat.shape[0]): #to exclude itself in minimum distance calcuration in family_distance
            self.distance_mat[i][i] = 9999
    
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
