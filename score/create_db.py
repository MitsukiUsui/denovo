#!/usr/bin/env python3

import sys
import re
import pandas as pd
from io import StringIO
from ScoreDbController import ScoreDbController

class ProdigalScore:
    """
    Class to handle prodigal supplemenatl score information
    """
    
    def __init__(self, supFilepath):
        self.supFilepath = supFilepath
        self.df_lst, self.seqhdr_lst = self.read_sup(supFilepath)
        assert len(self.df_lst) == len(self.seqhdr_lst)
        for i in range(len(self.df_lst)):
            self.df_lst[i] = self.assign_chunk(self.df_lst[i])
    
    def read_sup(self, supFilepath):
        """
        return list of DataFrame divided by seqhdr and corresponding seqhdr_lst
        """
        
        def todf(line_lst):
            return pd.read_csv(StringIO("".join(line_lst)), delimiter='\t')
        
        def seqhdr(line):
            m = re.findall("seqhdr=\"(.*?)\"", line)
            assert len(m) == 1
            return m[0].split()[0]

        df_lst = []
        seqhdr_lst = []
        line_lst = [] #record lines for 1 dataframe
        with open(supFilepath, "r") as f:
            for line in f:
                if len(line) > 0:
                    if "# Sequence Data:" in line: #the first header came
                        if len(line_lst) > 0:
                            df_lst.append(todf(line_lst))
                            line_lst = []
                        seqhdr_lst.append(seqhdr(line))
                    elif "# Run Data:" in line: #the second header came
                        pass
                    else:
                        line_lst.append(line)
        df_lst.append(todf(line_lst))
        return df_lst, seqhdr_lst
    
    def assign_chunk(self, df):
        """
        assign chunk_id column
        1 chunk correspond to a series of orfs which share same stop codon
        """

        col = "End" if df["Std"].iloc[0] == "+" else "Beg"
        prv = df[col].iloc[0]

        stopId = 0
        stopId_lst = [] # list of stop codon id
        for _, row in df[["Beg", "End", "Std"]].iterrows():
            if row[col] == prv:
                stopId_lst.append(stopId)
            else:
                stopId += 1
                stopId_lst.append(stopId)
                #start new chank
                col = "End" if row["Std"] == "+" else "Beg"
                prv = row[col]
        df["chunk_id"] = stopId_lst
        return df

def main(supFilepath, dbFilepath):
    print("START: create database from {}".format(supFilepath))
    ps =  ProdigalScore(supFilepath)
    sdc = ScoreDbController(dbFilepath)
    for df, table in zip(ps.df_lst, ps.seqhdr_lst):
        sdc.create(table)
        df.to_sql(table, sdc.con, if_exists="append", index=False)
    print("DONE: create {} table in {}".format(len(ps.df_lst), dbFilepath))

if __name__=="__main__":
    target = sys.argv[1]
    strain = sys.argv[2]
    supFilepath="/data/mitsuki/data/denovo/{}/annotation/prodigal/sup/{}.sup".format(target, strain)
    dbFilepath=supFilepath.replace(".sup", ".sq3")
    main(supFilepath, dbFilepath)