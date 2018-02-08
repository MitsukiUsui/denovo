#!/usr/bin/env python3

import sys
import numpy as np
from ete3 import Tree

sys.path.append("../helper")
from myio import get_strain_lst, get_distance_mat

def main(target):
    strain_lst = get_strain_lst(target)

    # filter strain greedy
    distance_mat = get_distance_mat(target)
    thres = 0.5
    msk = np.ones(len(strain_lst)).astype(bool)
    for i in range(len(strain_lst)):
        if msk[i]:
            for j in range(i + 1, len(strain_lst)):
                if distance_mat[i][j] < thres:
                    msk[j] = False
    print("DONE: filter to {} from {} strains".format(msk.sum(), len(msk)))
    
    outFilepath = "../data/{}/strain.fil.lst".format(target)
    with open(outFilepath, "w") as f:
        for i in range(len(strain_lst)):
            if msk[i]:
                f.write("{}\n".format(strain_lst[i]))
    print("DONE: output {}".format(outFilepath))

if __name__=="__main__":
    target = sys.argv[1]
    main(target)