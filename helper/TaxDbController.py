from myutil.myutil import DbController

class TaxDbController(DbController):
    def __init__(self, dbFilepath):
        super().__init__(dbFilepath)

    def accession2taxid(self, accessionVersion):
        query = 'SELECT taxid FROM accession2taxid WHERE accession_version = "{}"'.format(accessionVersion)
        success = self.execute(query)
        if success:
            ret=self.cur.fetchone()
            if ret is not None:
                return ret["taxid"]
        return 0
    
if __name__=='__main__':
    tdc = TaxDbController("/data/mitsuki/data/refseq/accession2taxid/prot.accession2taxid.sq3")
    print(tdc.accession2taxid("WP_091645391.1"))
    print(tdc.accession2taxid("hoge"))
