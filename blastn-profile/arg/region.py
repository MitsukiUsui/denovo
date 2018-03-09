#!/usr/bin/env python3

import sys
sys.path.append("../helper")
from myio import get_strain_lst

target = sys.argv[1]
strain_lst = get_strain_lst(target)
for strain in strain_lst:
    print("{},{}".format(target, strain))

