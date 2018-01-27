#!/usr/bin/env python3
import sys
sys.path.append("/home/mitsuki/software")
from myutil.myutil import myrun

#for strand, phase in [("+", 1), ("+", 2), ("-", 0), ("-", 1), ("-", 2)]:
for strand, phase in [("+", 0)]:
    for simId in range(100):
        for dist in range(1,10):
            dist /= 10
            seqName = "sim{}{}_{}_{}".format(strand, phase, simId, dist)
            cmd = "./mafft.sh {}".format(seqName)
            myrun(cmd)
