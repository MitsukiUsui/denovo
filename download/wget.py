#!/usr/bin/env python3

import sys
import pandas as pd
sys.path.append("../helper")
from myutil.myutil import myrun

def main(catalogFilepath):    
    directory_lst = ["/data/mitsuki/data/refseq/dnaseq",
                        "/data/mitsuki/data/refseq/gff",
                        "/data/mitsuki/data/refseq/cds",
                        "/data/mitsuki/data/refseq/faa"]
    suffix_lst = ["_genomic.fna.gz",
                    "_genomic.gff.gz",
                    "_cds_from_genomic.fna.gz",
                    "_protein.faa.gz"]
    
    catalog_df=pd.read_csv(catalogFilepath, sep="\t")
    for _, row in catalog_df.iterrows():
        for directory, suffix in zip(directory_lst, suffix_lst):
            filename=row["ftp_path"].split("/")[-1] + suffix
            ftpFilepath="{}/{}".format(row["ftp_path"], filename)
            outFilepath="{}/{}".format(directory, filename)
            cmd = "wget -q -O {} {}".format(outFilepath, ftpFilepath)
#            print(cmd)
            success = myrun(cmd)
            if success:
                print("DONE: download to {}".format(outFilepath))
            if not(success):
                print("ERROR: wget failed for {}".format(ftpFilepath))    

if __name__ == "__main__":
    target=sys.argv[1]
    catalogFilepath="../data/{}/catalog.tsv".format(target)
    main(catalogFilepath)
