#!/usr/bin/env python3

"""
create cluster.tsv according to output of sonicparanoid
"""

import sys
import re
import pandas as pd
import numpy as np
from Bio import SeqIO

sys.path.append("../helper")
from myio import get_strain_lst

def find_no_ortho(seqFilepath, sonic_srs):
    # get list of orf_id first
    orfId_lst=[]
    for rec in SeqIO.parse(seqFilepath, "fasta"):
        orfId_lst.append(rec.id)
    orfId_set=set(orfId_lst)
    assert len(orfId_lst) == len(orfId_set) #orf_id given in refseq.py need to be unique

    # get list of orf_id which is in sonicparanoid output (= has orthlog famimly)
    hasOrtho_set=set()
    for orfIds in sonic_srs.dropna():
        for orfId in orfIds.split(','):
            hasOrtho_set.add(orfId)

    assert orfId_set >= hasOrtho_set
    noOrtho_set = set(orfId_lst) - hasOrtho_set
    return list(noOrtho_set)
 
def parse_sonic(sonicFilepath, strain_lst):
    """
    parse output of sonicparanoid and return sonic_df which has lineage, size, and strain columns
    """
    
    def trimmer(txt):
        """
        trim score notation of given by sonic paranoid
            e.g.) ${orf_id_a},${orf_id_b}:${orf_id_score} -> ${orf_id_a},${orf_id_b}
        """
        ret_lst = []
        ptn = r"([^,]*):[.0-9]*"
        for block in txt.split(','):
            ret_lst.append(re.sub(ptn, r"\1", block))
        return ','.join(ret_lst)
    
    sonic_df=pd.read_csv(sonicFilepath, sep="\t")
    sonic_df = sonic_df.rename(columns = {"sp_in_grp": "lineage", "group_size": "size"})
    for strain in strain_lst:
        ret_lst=[]
        for txt in sonic_df[strain]:
            if txt == "*": #sonicparanoid assign * if no gene belongs to the family
                ret_lst.append(np.nan)
            else:
                ret_lst.append(trimmer(txt))
        sonic_df[strain] = ret_lst
    return sonic_df[["lineage", "size"] + strain_lst]

def main(sonicDirec, outFilepath, strain_lst):
    sonicFilepath="{}/out/multi_species/multispecies_clusters_all.tsv".format(sonicDirec)
    sonic_df = parse_sonic(sonicFilepath, strain_lst)
    
    dct_lst=[]
    dct_lst += sonic_df.to_dict("record")
    print("DONE: add {} ortholog family".format(len(dct_lst)))

    print("START: add non-orthlog family")
    for strain in strain_lst:
        seqFilepath="{}/in/{}".format(sonicDirec, strain)
        noOrtho_lst = find_no_ortho(seqFilepath, sonic_df[strain])
        for noOrtho in noOrtho_lst:
            dct = {strain: noOrtho,
                    "lineage": 1,
                    "size": 1}
            dct_lst.append(dct)
        print("\tDONE: add {} singletons from {}".format(len(noOrtho_lst), strain))
        
    out_df=pd.DataFrame(dct_lst)
    out_df["family"]=["family{}".format(i) for i in out_df.index]
    out_df=out_df[["family", "lineage", "size"] + strain_lst]
    out_df.index.name="ClusterID"
    out_df.to_csv(outFilepath, sep="\t")
    print("DONE: output to {}".format(outFilepath))

if __name__=="__main__":
    target=sys.argv[1]
    annotType=sys.argv[2]
    sonicDirec="/data/mitsuki/data/denovo/{}/annotation/{}/sonic".format(target, annotType)
    outFilepath="../data/{}/cluster.tsv".format(target)
    strain_lst = get_strain_lst(target)
    main(sonicDirec, outFilepath, strain_lst)
