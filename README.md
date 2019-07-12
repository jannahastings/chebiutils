# ChebiDbUtils package

This is a simple Python implementation of some functionality for locally caching and searching [ChEBI](http://www.ebi.ac.uk/chebi) data. It works from the ChEBI OBO download file (nightly build) as input and builds a rapid-access Python dictionary of entity data as well as a searchable Whoosh index of commonly searched fields. 

To setup the local cache of ChEBI data, run:

~~~~ 

from chebidblite import setupdb

setupdb.prepareCacheAndIndex() 

~~~~

The cache and index will be built into the current working directory. This only needs to be executed once, but can be re-executed if a newer version of ChEBI is required.

To perform searches (for example), use: 

~~~~ 

from chebidblite import searcher 
     
chebisearcher = searcher.ChebiSearcher()

res = chebisearcher.findAllChildrenOf("CHEBI:15377")

ids = [r.chebi_id for r in res]

names = [r.chebi_name for r in res]

~~~~

Just to access the stored data by chebi id (for example), use: 

~~~~ 

from chebidblite import dblite

dbchebi = dblite.ChebiDbLite()
dbchebi.initialize()

water = dbchebi.getEntity("CHEBI:15377")
print(water.chebi_name)
print(water.is_a)
print([dbchebi.getEntity(e).chebi_name for e in water.is_a])

~~~~

This library is in a very preliminary stage of development! 


