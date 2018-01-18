import pandas as pd
from io import StringIO

def read_gff(gffFilepath, additional_lst=None):
    column_lst=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    def parse_attribute(att):
        dct={}
        att_lst=att.split(";")
        for e in att_lst:
            if len(e.split("=")) == 2:
                k,v = e.split("=")
                dct[k] = v
        return dct

    with open(gffFilepath, "r") as f:
        line_lst=[]
        for line in f:
            if line[0]!="#":
                line_lst.append(line)
    gff_df=pd.read_csv(StringIO("".join(line_lst)), delimiter='\t', header=None, names=column_lst)
    if additional_lst is not None:
        att_df=pd.DataFrame([parse_attribute(att) for att in gff_df["attribute"]])
        try:
            gff_df=pd.concat([gff_df, att_df[additional_lst]], axis=1)
        except KeyError:
            print("ERROR: {} does not contain information on {}".format(gffFilepath, additional_lst))
    return gff_df

def write_gff(outFilepath, gff_df):
    column_lst=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    gff_df=gff_df[column_lst]
    gff_df.to_csv(outFilepath, index=False, header=False, sep="\t")
