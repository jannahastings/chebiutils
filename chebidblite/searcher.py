import itertools
import os

import whoosh.index as index
from whoosh.qparser import QueryParser

from .dblite import ChebiDbLite
from .dblite import ChebiEntity

class ChebiSearcher:
    """A search interface wrapping queries on a Whoosh index """
    
    IS_A = "is_a:"
    HAS_STRUC = "has_struc:True"
    HAS_ROLE = "has_role:"
    LEAF_NODE = "leaf_node:True"
    AND = " and "
    OR = " or "
    
    CHEM_ROOT = "CHEBI:24431"
    ROLE_ROOT = "CHEBI:50906"
    GROUP_ROOT = "CHEBI:24433"
    MOLENT_ROOT = "CHEBI:23367"
    
    def __init__(self,indexdir = "indexdir"):
        self.cacheDir = os.getenv('CHEBIDBLITECACHE', '~')
        self.indexdir = indexdir
        self.ix = index.open_dir(self.cacheDir+self.indexdir)
        self.parser = QueryParser("chebi_name", self.ix.schema)
        self.db = ChebiDbLite()
        self.db.initialize()
        
    def _processSearchResultsSingle(self,results):
        if len(results) > 0:
            chebiId = results[0]["chebi_id"]
            if chebiId in self.db.data_dict.keys():
                return (self.db.data_dict[chebiId])
            else: 
                print("CHEBI ID: ",chebiId,"NOT FOUND")
        else:
            return (None)
        
    def _processSearchResultsList(self,results):
        results_ids = set()
        for hit in results:
            results_ids.add(hit["chebi_id"])
        if len(results_ids)>0:
            return (set([self.db.data_dict[x] for x in results_ids]))
        else:
            return (None)
        
    def findChebiIdByName(self,name):
        with self.ix.searcher() as searcher:
            query = self.parser.parse(name)
            results = searcher.search(query,limit=1)
            return (self._processSearchResultsSingle(results))
        
    def findAllChildrenOf(self,chebiId):
        with self.ix.searcher() as searcher:
            query = self.parser.parse(self.IS_A+chebiId)
            results = searcher.search(query,limit=None)
            return (self._processSearchResultsList(results))
                
    def findAllChildrenWithStructures(self,chebiId):
        with self.ix.searcher() as searcher:
            query = self.parser.parse(self.IS_A+chebiId + self.AND + self.HAS_STRUC)
            results = searcher.search(query,limit=None)
            return (self._processSearchResultsList(results))
     
    def findAllLeafChildrenWithStructures(self,chebiId):
        with self.ix.searcher() as searcher:
            query = self.parser.parse(self.IS_A+chebiId + self.AND + self.HAS_STRUC + self.AND + self.LEAF_NODE)
            results = searcher.search(query,limit=None)
            return (self._processSearchResultsList(results))
                
    def findAllChildrenWithRole(self,chebiId,roleId):
        with self.ix.searcher() as searcher:
            query = self.parser.parse(self.IS_A+chebiId + self.AND + self.HAS_ROLE + roleId)
            results = searcher.search(query,limit=None)
            return (self._processSearchResultsList(results))

    def findAllWithRole(self,roleId):
        return (self.findAllChildrenWithRole(self.CHEM_ROOT,roleId))
        
    def findAllByQueryString(self, queryString):
        with self.ix.searcher() as searcher:
            query = self.parser.parse(queryString)
            results = searcher.search(query,limit=None)
            return (self._processSearchResultsList(results))
