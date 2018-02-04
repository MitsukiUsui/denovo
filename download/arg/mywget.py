#!/usr/bin/env python3

import sys
import pandas as pd

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
            outFilepath="{}/{}".format(directory, filename.replace(".gz", ""))
            print("{},{}".format(ftpFilepath, outFilepath))

if __name__ == "__main__":
    target=sys.argv[1]
    catalogFilepath="../data/{}/catalog.tsv".format(target)
    main(catalogFilepath)
