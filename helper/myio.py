import numpy as np
from ete3 import Tree

direc="/home/mitsuki/altorf/denovo/data"

def get_strain_lst(target, filtered = False):
    if filtered:
        fp = "{}/{}/strain.fil.lst".format(direc, target)
    else:
        fp = "{}/{}/strain.lst".format(direc, target)
    return [s.strip() for s in open(fp, 'r').readlines()]

def get_distance_mat(target):
    fp = "{}/{}/distance.npy".format(direc, target)
    try:
        distance_mat = np.load(fp)
    except FileNotFoundError:
        distance_mat = calc_distance_mat(target)
        np.sace(fp, distance_mat)
    return distance_mat

def calc_distance_mat(target):
    fp = "{}/{}/cluster.phb".format(direc, target)
    tree = Tree(fp)

    strain_lst = get_strain_lst(target)
    size = len(strain_lst)
    distance_mat = -np.ones((size, size))
    for i in range(size):
        for j in range(i, size):
            distance_mat[i][j] = tree.get_distance(strain_lst[i], strain_lst[j])
    return distance_mat

if __name__=="__main__":
    print(get_strain_lst("synechococcus"))
