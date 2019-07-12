import itertools
import os
from whoosh.index import create_in
from whoosh.fields import *
from .dblite import ChebiDbLite
from .dblite import ChebiEntity

class ChebiIndexer:
    """This class creates a simple Whoosh index of ChEBI data"""
    
    schema = Schema(chebi_id=ID(stored=True),
                chebi_name=TEXT(phrase=True),
                definition=TEXT(),
                has_struc=BOOLEAN,
                struc=TEXT(phrase=True),
                leaf_node=BOOLEAN,
                is_a = TEXT,
                has_role = TEXT 
                )
    
    def __init__(self,indexdir="indexdir"):
        self.cacheDir = os.getenv('CHEBIDBLITECACHE', '~')
        self.indexdir = indexdir
        self.ix = None
        self.db = ChebiDbLite()
        self.db.initialize()
        self.debug = True
        

    def buildIndex(self):
        # Need to do this only if needed
        os.makedirs(self.cacheDir+self.indexdir,exist_ok=True)

        # Need to sanity check the state (maps already built successfully etc)?
        if not self.db.initialized:
            print("ERROR: Data store is not initialized. Exiting.")
            return()

        # Build the index    TODO check if locked, work around
        self.ix = create_in(self.cacheDir+self.indexdir, self.schema)
        writer = self.ix.writer(limitmb=256)
        try: 
            if (self.debug): 
                print("Writing to index...")
            # use ChebiEntity to store data. todo expand indexed fields for searching
            i = 1
            for entity in self.db.data_dict.values():
                i = i+1
                if self.debug:
                    if i % 10000 == 0:
                        print("Writing entity no. ",i)
                leaf_node = False
                if not entity.children or len(entity.children) == 0:
                    leaf_node = True
                writer.add_document(chebi_id=entity.chebi_id,
                                    chebi_name=entity.chebi_name,
                                    definition=entity.definition,
                                    has_struc=entity.has_struc,
                                    struc=entity.smiles,
                                    leaf_node = leaf_node,
                                    is_a = " ".join(self.db.ancestor_map[entity.chebi_id]),
                                    has_role = " ".join(self.db.role_map[entity.chebi_id]))
            if (self.debug): 
                print("Committing (saving) index...")
            writer.commit()
        except Exception:
            exc_type, value, traceback = sys.exc_info()
            print("Something went wrong:",exc_type,value,traceback)
            writer.cancel()

    

    def getIndex(self):
        return self.ix










