#!/usr/bin/env python3
import sys
sys.path.append("/home/mitsuki/software")
from myutil.myutil import myrun

for strand, phase in (("+", 1), ("+", 2), ("-", 0), ("-", 1), ("-", 2)):
    for dist in range(1, 10):
        dist /= 10
        seqName = "sim{}{}_0_{}".format(strand, phase, dist)
        for trainName in ("train2000", "train10-4", "train10-5", "train10-6"):
            cmd = "./caller.sh {} {}".format(seqName, trainName)
            myrun(cmd)
#            print(cmd)
