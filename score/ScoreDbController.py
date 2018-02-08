import sys
sys.path.append("../helper")

from myutil.myutil import DbController

class ScoreDbController(DbController):
    def __init__(self, dbFilepath):
        super().__init__(dbFilepath)
        
    def create(self, table):
        if table in self.table_lst:
            self.clear_row(table)
            print("DEBUG: clear row from {}".format(table), file=sys.stderr)
        else:
            query = ('CREATE TABLE "{}" ('.format(table)\
                     +'"Beg" INTEGER, "End" INTEGER, "Std" TEXT, "Total" REAL, "CodPot" REAL, "StrtSc" REAL, "Codon" TEXT, "RBSMot" TEXT, "Spacer" TEXT,'\
                     +'"RBSScr" REAL, "UpsScr" REAL,"TypeScr" REAL,"GCCont" REAL,"chunk_id" INTEGER, PRIMARY KEY (Beg, End));')
            success = self.execute(query)
            if success:
                self.table_lst.append(table)
                print("DONE: create table {}".format(table), file=sys.stderr)
            else:
                print("ERROR: failed to create table {}".format(table), file=sys.stderr)
            
    def info(self, seqhdr, begin, end):
        assert seqhdr in self.table_lst
        
        query = ('SELECT * FROM "{0}" WHERE Beg = {1} AND End = {2}'.format(seqhdr, begin, end))
        success = self.execute(query)
        if success:
            ret = self.cur.fetchone()
            if ret is None:
#                print("ERROR: {}-[{},{}] does not exist".format(seqhdr, begin, end), file=sys.stderr)
                return dict()
            else:
                return ret
        else:
            return dict()
