#!/usr/bin/env python3

import sys

target=sys.argv[1]
strainFilepath="../data/{}/strain.lst".format(target)
strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
for strain in strain_lst:
    print("{},{}".format(target, strain))
