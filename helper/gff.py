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
        assert att_df.shape[0] == gff_df.shape[0]
        for col in additional_lst:
            try:
                gff_df[col] = att_df[col]
            except KeyError:
                print("ERROR: {} does not contain information on {}".format(gffFilepath, col))
    return gff_df

def write_gff(gff_df, outFilepath):
    column_lst=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    gff_df=gff_df[column_lst]
    gff_df.to_csv(outFilepath, index=False, header=False, sep="\t")

if __name__ == "__main__":
    gffFilepath = "/data/mitsuki/data/refseq/gff/GCF_000011345.1_ASM1134v1_genomic.gff"
    gff_df = read_gff(gffFilepath, ['ID', 'Parent', 'locus_tag', 'protein_id', 'pseudo'])
    print(gff_df.columns)
