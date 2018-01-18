#!/usr/bin/env python3

import sys
sys.path.append("/home/mitsuki/altorf/mbgd/helper")
from myutil.myutil import myrun

target=sys.argv[1]
strainFilepath="../../data/{}/strain.lst".format(target)
strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
for strain in strain_lst:
    cmd = "qsub caller.sh {} {}".format(target, strain)
    myrun(cmd)
