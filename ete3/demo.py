#!/usr/bin/env python
from ete3 import PhyloTree

fp="demo.nwk"
ofp="demo.png"
t = PhyloTree(fp)
t.render(ofp, w=2000, dpi=500)