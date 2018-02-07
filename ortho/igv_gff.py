#!/usr/bin/env python3

import sys
import os
sys.path.append("../helper")
from gff import read_gff, write_gff
from myutil.myutil import myrun

def format_attribute(att_str, delCol_lst = None, addCol_dct = None):
    """
    delCol_lst = [delCol1, delCol2, ...]
    addCol_dct = {"newCol1": newCol1Val, "newCol2": newCol2Val, ...}
    """
    
    def str2dct(att_str):
        dct={}
        for e in att_str.split(";"):
            if len(e.split("=")) == 2:
                k,v = e.split("=")
                dct[k] = v
        return dct
    
    def dct2str(att_dct):
        ret_lst = []
        for key,val in att_dct.items():
            ret_lst.append("{}={}".format(key, val))
        return ';'.join(ret_lst)
    
    att_dct = str2dct(att_str)
    
    if delCol_lst is not None:
        for col in delCol_lst:
            try:
                del att_dct[col]
            except KeyError:
                pass
    
    if addCol_dct is not None:
        att_dct.update(addCol_dct)        
    
    return dct2str(att_dct)


def refseq(inFilepath, outFilepath):
    """
    delete some information in attribute, and renew ID column with unique id
    """
    
    refseq_df=read_gff(inFilepath, ["family", "orf_id"])
    att_lst = []
    for _, row in refseq_df.iterrows():
        delCol_lst=["Parent", "ID", "gene"]
        addCol_dct = {"ID": "{}({})".format(row["family"],row["orf_id"])}
        att_lst.append(format_attribute(row["attribute"], delCol_lst = delCol_lst, addCol_dct = addCol_dct))
    refseq_df["attribute"]=att_lst
    write_gff(outFilepath, refseq_df)
    
def prodigal(inFilepath, outFilepath):
    """
    just create symbolic link so far
    """
    
    cmd = "ln -s {} {}".format(inFilepath, outFilepath)
    myrun(cmd)
    
def main(target, strain_lst):
    outDirec="/data/mitsuki/data/denovo/{}/igv".format(target)
    os.makedirs(outDirec, exist_ok=True)
    
    for strain in strain_lst:
        inRefFilepath="/data/mitsuki/data/denovo/{}/annotation/refseq/gff/{}.gff".format(target, strain)
        outRefFilepath="{}/{}_refseq.gff".format(outDirec, strain)
        refseq(inRefFilepath, outRefFilepath)

        inProFilepath="/data/mitsuki/data/denovo/{}/annotation/prodigal/gff/{}.gff".format(target, strain)
        outProFilepath="{}/{}_prodigal.gff".format(outDirec, strain)
        prodigal(inProFilepath, outProFilepath)
    print("DONE: output .gffs to {}".format(outDirec))

if __name__=="__main__":
    target=sys.argv[1]
    strainFilepath="../data/{}/strain.lst".format(target)
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    main(target, strain_lst)

    
    