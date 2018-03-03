import pandas as pd
import numpy as np
from ete3 import Tree

direc="/home/mitsuki/altorf/denovo/data"

def get_strain_lst(target, full=False):
    catalogFilepath = "{}/{}/catalog.tsv".format(direc, target)
    catalog_df = pd.read_csv(catalogFilepath, sep="\t")
    if full:
        return list(catalog_df["genome_id"])
    else:
        return list(catalog_df[catalog_df["represent"]==1]["genome_id"])

def get_cluster_df(target):
    fp = "{}/{}/cluster.tsv".format(direc, target)
    df = pd.read_csv(fp, delimiter='\t', dtype=object)
    df["lineage"] = df["lineage"].astype(int)
    df["size"] = df["size"].astype(int)
    return df

def get_distance_mat(target, full=False):
    fp = "{}/{}/distance.npy".format(direc, target)
    try:
        distance_mat = np.load(fp)
    except FileNotFoundError:
        distance_mat = calc_distance_mat(target)
        np.save(fp, distance_mat)

    if full:
        return distance_mat
    else:
        catalogFilepath = "{}/{}/catalog.tsv".format(direc, target)
        catalog_df = pd.read_csv(catalogFilepath, sep="\t")
        msk = np.array(catalog_df["represent"] == 1)
        return distance_mat[msk, :][:, msk]
        

def calc_distance_mat(target):
    fp = "{}/{}/cluster.phb".format(direc, target)
    tree = Tree(fp)

    strain_lst = get_strain_lst(target, full=True)
    size = len(strain_lst)
    distance_mat = -np.ones((size, size))
    for i in range(size):
        for j in range(i + 1, size):
            distance_mat[i][j] = tree.get_distance(strain_lst[i], strain_lst[j])
            distance_mat[j][i] = distance_mat[i][j]
    return distance_mat

if __name__=="__main__":
    print(get_strain_lst("synechococcaceae"))
