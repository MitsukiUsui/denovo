from ete3 import NCBITaxa

class NCBIController:
    def __init__(self):
        self.ncbi = NCBITaxa()
        
    def translate(self, taxid):
        """
        :ret scientific name
        """
        return self.ncbi.get_taxid_translator([taxid])[taxid]
    
    def get_lineage(self, taxid, rank_lst = None):        
        if rank_lst is None:
            rank_lst = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

        dct={}
        try:
            for taxidLineage, rank in self.ncbi.get_rank(self.ncbi.get_lineage(taxid)).items():
                if rank in rank_lst:
                    dct[rank] = taxidLineage
                    dct[rank+"_s"] = self.translate(taxidLineage)
            return dct
        except (KeyError, ValueError):
#            print("ERROR: unknown taxid = {}".format(taxid))
            return dict()
        
    def get_descendant(self, taxid, rank):
        ret = []
        children = self.ncbi.get_descendant_taxa(taxid, rank_limit=rank)
        for k, v in self.ncbi.get_rank(children).items():
            if v == rank:
                ret.append(k)
        return ret
    
if __name__=="__main__":
    nc = NCBIController()
    print(nc.get_descendant(85025, rank="genus"))
