#!/usr/bin/env python3

import sys
sys.path.append("../helper")

from myio import get_strain_lst

target=sys.argv[1]
strain_lst = get_strain_lst(target, full=False)
for strain in strain_lst:
    filepath = "/data/mitsuki/data/denovo/{}/dnaseq/{}.dnaseq".format(target, strain)
    query = "./query/{}/{}.fna".format(target, strain)
    result = "./result/{}/{}.tsv".format(target, strain)
    print("{},{},{},{}".format(strain, filepath, query, result))
