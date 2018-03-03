#!/usr/bin/env python

import sys
from ete3 import PhyloTree

target=sys.argv[1]
nwkFilepath="/home/mitsuki/altorf/denovo/data/{}/cluster.format.phb".format(target)
outFilepath="{}.png".format(target)
t = PhyloTree(nwkFilepath)
t.render(outFilepath, w=2000, dpi=500)
