#!/usr/bin/env python3

import sys
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import Counter

sys.path.append("../helper")
from gff import read_gff, write_gff

def get_fp(target, basename, genomeId):
    """
    organize various filepath information into infp & outfp respectively
    basename: name defined in ftp_path column
    """
    
    infp={}
    direc="/data/mitsuki/data/refseq"
    name=basename
    infp["dnaseq"]="{}/dnaseq/{}_genomic.fna".format(direc, name)
    infp["fna"]="{}/cds/{}_cds_from_genomic.fna".format(direc, name)
    infp["faa"]="{}/faa/{}_protein.faa".format(direc, name)
    infp["gff"]="{}/gff/{}_genomic.gff".format(direc, name)

    outfp={}
    direc="/data/mitsuki/data/denovo/{}/annotation/refseq".format(target)
    name=genomeId
    outfp["dnaseq"]="/data/mitsuki/data/denovo/{}/dnaseq/{}.dnaseq".format(target, name) # only dnaseq is stored upper directory
    outfp["fna"]="{}/fna/{}.fna".format(direc, name)
    outfp["faa"]="{}/faa/{}.faa".format(direc, name)
    outfp["gff"]="{}/gff/{}.gff".format(direc, name)
        
    return infp, outfp

def edit_dnaseq(inFilepath, outFilepath, genomeId):
    """
    add genome_id to head of seqname
    """
    with open(outFilepath, "w") as f:
        for record in SeqIO.parse(inFilepath, "fasta"):
            record.id = "{}-{}".format(genomeId ,record.id)  #add  genome_id to head
            SeqIO.write(record,  f, "fasta")

def edit_gff(gff_df, outFilepath):
    """
    add genome_id to head of seqname, add orf_id attribute
    """
    assert "orf_id" in gff_df.columns
    
    gff_df = gff_df[gff_df["feature"]=="CDS"].copy()  #filter only CDS
    gff_df["seqname"] = ["{}-{}".format(genomeId, seqname) for seqname in gff_df["seqname"]]

    att_lst = []
    for _, row in gff_df.iterrows():
        att = "{};orf_id={}".format(row["attribute"], row["orf_id"])
        att_lst.append(att)
    gff_df["attribute"]=att_lst
    write_gff(gff_df, outFilepath)
       
def parse_description(description):
    """
    parse description straing written in .fna
    """
    element_lst = re.findall('\[(.*?)\]', description) # ? for shortest match
    dct = {}
    for e in element_lst:
        e = e.split("=")
        if len(e)==2:
            dct[e[0]]=e[1]
    return dct

def set_orfId(gff_df):
    """
    format of orf_id is ${genome_id}-${cds_id}_${cds_id_count}
    cds_id_count is required only when cds_id is not unique
    """
    assert "ID" in gff_df.columns
    
    counter = Counter(gff_df[gff_df["feature"]=="CDS"]["ID"])
    duplication_lst=[key for key, val in counter.items() if val > 1]
#    print("DEBUG: found {} duplication for cds_id".format(len(duplication_lst), file=sys.stderr))
    
    #assign orfId for every row based on "ID" and counter
    orfId_lst=[] 
    for _, row in gff_df.iterrows():
        if row["feature"] != "CDS":
            orfId_lst.append(np.nan)
        else:
            cdsId = row["ID"]
            if cdsId in duplication_lst:
                orfId = "{}-{}_{}".format(genomeId, cdsId, counter[cdsId])
                counter[cdsId]-=1
            else:
                orfId = "{}-{}".format(genomeId, cdsId)
            orfId_lst.append(orfId)
    gff_df["orf_id"] = orfId_lst
    return gff_df

def get_lookup_df(gff_df): 
    for col in ["ID","Parent", "locus_tag", "protein_id"]:
        assert col in gff_df.columns
    if "pseudo" not in gff_df.columns:
        gff_df["pseudo"] = False
    
    gene2loc={} # given gene_id, return locus_tag
    dct_lst=[]
    for _, row in gff_df.iterrows():
        if row["feature"] == "gene":
            gene2loc[row["ID"]] = row["locus_tag"]
        elif row["feature"]=="CDS":
            dct = {}
            dct["orf_id"] = row["orf_id"]
            dct["protein_id"] = row["protein_id"]

            if row["Parent"] in gene2loc.keys(): #assume Parent gene record is already processed
                dct["locus_tag"] = gene2loc[row["Parent"]]
            elif row["pseudo"] == "true":
                pass
            else:
                print("WARN: {} lack protein information for some reason.".format(dct["orf_id"]), file = sys.stderr)
            dct_lst.append(dct)
    lookup_df=pd.DataFrame(dct_lst)
    lookup_df=lookup_df[["orf_id", "locus_tag", "protein_id"]]
    return lookup_df
        
def main(target, basename, genomeId):
    
    #organize input & output filepath
    infp, outfp = get_fp(target, basename, genomeId)
    print("IN:")
    for category in ("dnaseq", "fna", "faa", "gff"):
        print("\t{}".format(infp[category]))
    print("OUT:")
    for category in ("dnaseq", "fna", "faa", "gff"):
        print("\t{}".format(outfp[category]))
    print()
        
    
    edit_dnaseq(infp["dnaseq"], outfp["dnaseq"], genomeId)
    print("DONE: output to {}".format(outfp["dnaseq"]))
    
    gff_df = read_gff(infp["gff"], ["ID", "Parent", "locus_tag", "protein_id", "pseudo"])
    gff_df = set_orfId(gff_df)
    lookup_df=get_lookup_df(gff_df)
    
    print("START: load FASTA")
    loc2rec={} #locus_tag to dna record
    for rec in SeqIO.parse(infp["fna"], "fasta"):
        description = parse_description(rec.description)
        loc2rec[description["locus_tag"]] = rec
    print("\tDONE: load {} dna seq".format(len(loc2rec)))
    pro2rec={} #protein_id to protein record
    for rec in SeqIO.parse(infp["faa"], "fasta"):
        pro2rec[rec.id] = rec
    print("\tDONE: load {} protein seq".format(len(pro2rec)))
    
    print("DONE: output FASTA for {}/{} cds".format(lookup_df.dropna().shape[0], lookup_df.shape[0]))
    with open(outfp["fna"], "w") as fna, open(outfp["faa"], "w") as faa:
        for _, row in lookup_df.dropna().iterrows():
            dnaRec=loc2rec[row["locus_tag"]]
            desc = parse_description(dnaRec.description)
            if ("protein_id" in desc.keys()) and desc["protein_id"] != row["protein_id"]:  #check the validness of lookup_df
                print("ERROR: lookup table has incorrect matching for {}".format(row["orf_id"]), file=sys.stderr)
                exit(1)
            dnaRec.id = row["orf_id"]
            SeqIO.write(dnaRec, fna, "fasta")

            proRec = pro2rec[row["protein_id"]]
            proRec.id = row["orf_id"]
            SeqIO.write(proRec, faa, "fasta")
        print("\tDONE: output {}".format(outfp["fna"]))
        print("\tDONE: output {}".format(outfp["faa"]))

    edit_gff(gff_df, outfp["gff"])
    print("DONE: output to {}".format(outfp["gff"]))

if __name__=="__main__":
    target = sys.argv[1]
    basename = sys.argv[2]
    genomeId = sys.argv[3]
    main(target, basename, genomeId)
