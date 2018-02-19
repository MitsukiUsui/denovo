#!/usr/bin/env python

import sys
import pandas as pd
from ete3 import PhyloTree, Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace

sys.path.append("/home/mitsuki/altorf/denovo/helper")
from myio import *

def get_strain(cluster_df, strain_lst, family):
    return list(cluster_df[cluster_df["family"]==family][strain_lst].dropna(axis=1).columns)

def layout(node):
    global f1_lst
    global f2_lst

    if node.is_leaf():
        N = AttrFace("name", fsize=14, fgcolor="black") #  can be "sci_name"
        faces.add_face_to_node(N, node, 0)
        
        if node.name in f1_lst:
            C = CircleFace(radius=5, color="darkred", style="sphere")
            faces.add_face_to_node(C, node, 1, position="aligned")
        if node.name in f2_lst:
            C = CircleFace(radius=5, color="RoyalBlue", style="sphere")
            faces.add_face_to_node(C, node, 2, position="aligned")


target=sys.argv[1]
family1=sys.argv[2]
family2=sys.argv[3]
outFilepath="{}_{}_{}.png".format(target, family1, family2)
#outFilepath=sys.argv[4]

dataDirec="/home/mitsuki/altorf/denovo/data/{}".format(target)
nwkFilepath="{}/cluster.phb".format(dataDirec)
strain_lst = get_strain_lst(target)
cluster_df = get_cluster_df(target)

f1_lst=get_strain(cluster_df, strain_lst, family1)
f2_lst=get_strain(cluster_df, strain_lst, family2)

t = PhyloTree(nwkFilepath)
ts = TreeStyle()
ts.layout_fn = layout
ts.show_leaf_name = False
#t.show(tree_style=ts)
t.render(outFilepath, tree_style=ts, w=183, units="mm")
