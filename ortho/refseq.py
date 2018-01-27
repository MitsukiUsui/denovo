#!/usr/bin/env python3

import sys
import re
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

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
            record.id = "{}:{}".format(genomeId ,record.id)  #add  genome_id to head
            SeqIO.write(record,  f, "fasta")

def edit_gff(inFilepath, outFilepath, genomeId):
    """
    add genome_id to head of seqname, add orf_id attribute
    """
    gff_df=read_gff(inFilepath, ["ID", "protein_id"])
    gff_df = gff_df[gff_df["feature"]=="CDS"]  #filter only CDS
    gff_df["seqname"] = ["{}:{}".format(genomeId, seqname) for seqname in gff_df["seqname"]]
    gff_df["orf_id"] = ["{}:{}".format(genomeId, cds_id) for cds_id in gff_df["ID"]]

    att_lst = []
    for _, row in gff_df.iterrows():
        att = "{};orf_id={}".format(row["attribute"], row["orf_id"])
        att_lst.append(att)
    gff_df["attribute"]=att_lst
    write_gff(outFilepath ,gff_df)

def get_lookup_df(gffFilepath, fnaFilepath, faaFilepath):
    """
    connect 3 file and create lookup table
    orf_id, fna_id, protein_id is used as respective identifier.
    currently, faaFilepath is not used at all, but we need to check whether protein_id is valid or not before returning.
    """
    
    def parse_description(description):
        """
        parse description straing written in .fna
        """
        
        element_lst = re.findall('\[(.*?)\]', record.description) # ? for shortest match
        dct = {}
        for e in element_lst:
            e = e.split("=")
            if len(e)==2:
                dct[e[0]]=e[1]
        return dct
    
    print("START: connect 3 information written in gff, fna, and faa")
    
    gff_df=read_gff(gffFilepath, ["orf_id", "protein_id"])
    lookup_df1=gff_df[["orf_id", "protein_id"]]
    lookup_df1 = lookup_df1.rename(columns = {"protein_id": "protein_id1"})
    lookup_df1 = lookup_df1.drop_duplicates(subset=['orf_id'])
    print("\tDONE: find {} CDSs in .gff. remove {} duplication".format(lookup_df1.shape[0] , gff_df.shape[0] - lookup_df1.shape[0]))
    
    dct_lst=[]
    for record in SeqIO.parse(fnaFilepath, "fasta"):
        dct={}
        dct["fna_id"]=record.id
        dct["orf_id"]="{}:cds{}".format(genomeId, int(dct["fna_id"].split("_")[-1]) - 1)
        description = parse_description(record.description)
        if "protein_id" in description.keys():
            dct["protein_id2"] = description["protein_id"]
        dct_lst.append(dct)
    lookup_df2=pd.DataFrame(dct_lst)
    print("\tDONE: find {} dna seqs in .fna".format(lookup_df2.shape[0]))

    assert lookup_df1.shape[0] == lookup_df2.shape[0]
    lookup_df=pd.merge(lookup_df1, lookup_df2, on="orf_id", how="inner")
    assert lookup_df.shape[0] == lookup_df1.shape[0]
    pro_lst=[]
    for pro1, pro2 in zip(lookup_df["protein_id1"], lookup_df["protein_id2"]):
        if pro1 == pro2:
            pro_lst.append(pro1)
        else:
            pro_lst.append(None)
    lookup_df["protein_id"]=pro_lst
    lookup_df=lookup_df[["orf_id", "fna_id", "protein_id"]]
    lookup_df=lookup_df[lookup_df.isnull().sum(axis=1) == 0]
    print("\tDONE: find {} valid commbination".format(lookup_df.shape[0]))
    return lookup_df
    
def get_lookup_dct(lookup_df):
    fna2orf = {}
    pro2orf = defaultdict(list)
    for _, row in lookup_df.iterrows():
        fna2orf[row["fna_id"]] = row["orf_id"] # fna_id and orf_id one to one
        pro2orf[row["protein_id"]].append(row["orf_id"]) # multiple orf_id can have same fna_id
    return fna2orf, pro2orf
    
def edit_fna(inFilepath, outFilepath, fna2orf):
    with open(outFilepath, "w") as f:
        totalCount, outputCount = 0, 0
        for record in SeqIO.parse(inFilepath, "fasta"):
            totalCount+=1
            if record.id in fna2orf.keys():
                outputCount+=1
                record.id = fna2orf[record.id]
                record.description = ""
                SeqIO.write(record,  f, "fasta")
    print("DONE: output {}/{} dna seqs to {}".format(outputCount, totalCount, outFilepath))
    
def edit_faa(inFilepath, outFilepath, pro2orf):
    with open(outFilepath, "w") as f:
        totalCount, outputCount = 0, 0
        for record in SeqIO.parse(inFilepath, "fasta"):
            totalCount+=1
            for orfId in pro2orf[record.id]:
                outputCount += 1
                record.id = orfId
                record.description = "" #need to remove description because its too long for phylophlan (FastTree)
                SeqIO.write(record,  f, "fasta")
    print("DONE: output {}/{} protein seqs to {}".format(outputCount, totalCount, outFilepath))
    
def main(target, basename, genomeId):
    infp, outfp = get_fp(target, basename, genomeId)
    
    print("IN")
    for category in ("dnaseq", "fna", "faa", "gff"):
        print("\t{}".format(infp[category]))
    print("OUT")
    for category in ("dnaseq", "fna", "faa", "gff"):
        print("\t{}".format(outfp[category]))
    print()
        
    edit_dnaseq(infp["dnaseq"], outfp["dnaseq"], genomeId)
    print("DONE: output to {}".format(outfp["dnaseq"]))
    
    edit_gff(infp["gff"], outfp["gff"], genomeId)
    print("DONE: output to {}".format(outfp["gff"]))
    
    lookup_df=get_lookup_df(outfp["gff"], infp["fna"], infp["faa"])
    if lookup_df.shape[0] < 10:
        print("ERROR: seemed to fail to finding valid combination between gff, fna and faa.", file=sys.stderr)
        print("\tplease exclude {} ({}) from downstream analysis".format(genomeId, basename), file=sys.stderr)
        exit(1)
    else:
        fna2orf, pro2orf = get_lookup_dct(lookup_df)
        edit_fna(infp["fna"], outfp["fna"], fna2orf)
        edit_faa(infp["faa"], outfp["faa"], pro2orf)

if __name__=="__main__":
    target = sys.argv[1]
    basename = sys.argv[2]
    genomeId = sys.argv[3]
    main(target, basename, genomeId)
    