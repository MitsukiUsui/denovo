import random
import pandas as pd
import numpy as np
from ete3 import Tree
from Bio import SeqIO

direc="/home/mitsuki/altorf/denovo/data"

class DistanceManager:
    def __init__(self, target):
        distance_mat = get_distance_mat(target)
        for i in range(distance_mat.shape[0]):
            distance_mat[i][i] = 0
        self.distance_mat = distance_mat
        self.strain_lst = get_strain_lst(target)

    def distance(self, g1, g2):
        return self.distance_mat[self.strain_lst.index(g1)][self.strain_lst.index(g2)]

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

def get_orf2family(target):
    orf2family = {}
    cluster_df = get_cluster_df(target)
    strain_lst = get_strain_lst(target)
    for strain in strain_lst:
        for family, orfids in zip(cluster_df["family"], cluster_df[strain]):
            if isinstance(orfids, str):
                for orfid in orfids.split(','):
                    orf2family[orfid] = family
    return orf2family

def sample_record(target, family, ext, n=1):
    assert ext=="fna" or ext=="faa"
    fp = "/data/mitsuki/data/denovo/{0}/annotation/refseq/family/{2}/{1}.{2}".format(target, family, ext)
    rec_lst = []
    for rec in SeqIO.parse(fp, "fasta"):
        rec_lst.append(rec)
    return random.sample(rec_lst, min(len(rec_lst), n))

def get_hit_df(fp):
    columns_lst=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart","send", "evalue", "bitscore"]
    hit_df=pd.read_csv(fp, delimiter='\t', header=None, names=columns_lst)
    return hit_df

if __name__=="__main__":
    print(get_strain_lst("synechococcaceae"))
