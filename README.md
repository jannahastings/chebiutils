# ChebiDbUtils package

This is a simple implementation of functionality for locally caching and searching [ChEBI](http://www.ebi.ac.uk/chebi) data. It works from the ChEBI OBO download file (nightly build) as input and builds a rapid-access Python dictionary of entity data as well as a searchable Whoosh index of commonly searched fields. 

To setup the local cache of ChEBI data, run:

~~~~ from chebidblite import setupdb

setupdb.prepareCacheAndIndex() ~~~~

The cache and index will be built into the current working directory.

To perform searches or access cached data thereafter, use: 

~~~~ from chebidblite import searcher 
     
     chebisearcher = searcher.ChebiSearcher()

     res = chebisearcher.findAllChildrenOf("CHEBI:15377")

     ids = [r.chebi_id for r in res]

     names = [r.chebi_name for r in res]

~~~~

  
