# ChebiUtils 
## ChebiDbLite

This is a simple Python implementation of some functionality for locally caching and searching [ChEBI](http://www.ebi.ac.uk/chebi) data. It works from the ChEBI OBO download file (nightly build) as input and builds a rapid-access Python dictionary of entity data as well as a searchable Whoosh index of commonly searched fields. The main purpose of this library is to enable accessing ChEBI content and traversing common relational patterns *quickly*.

Beware: this library is in a very preliminary stage of development! 

To install, clone the repository to a directory on your local file system and install from there using pip:
~~~~
pip install -e /path/to/your/directory/chebiutils
~~~~

To setup the local cache of ChEBI data, run:

~~~~ 

from chebidblite import setupdb

setupdb.prepareCacheAndIndex() 

~~~~

The cache and index will be built into the directory specified by the environment variable CHEBIDBLITECACHE, which defaults to the following folder within your user home directory: '~/Library/Caches/chebidblite/'. This takes a few minutes. It only needs to be executed once, but can be re-executed whenever a newer version of ChEBI is required -- e.g., overnight.

Once the cache and index are prepared, they can be accessed from any other Python script as required without the overhead of rebuilding. 

To perform searches (for example), use: 

~~~~ 

from chebidblite import searcher 
     
chebisearcher = searcher.ChebiSearcher()

res = chebisearcher.findAllChildrenOf("CHEBI:15377")

ids = [r.chebi_id for r in res]

names = [r.chebi_name for r in res]

~~~~

Another example, looking for descendents of a particular class that have structures: 
~~~~
# get all leaf nodes that are descendents of 'tricarboxylic acid' and have structures
res = chebisearcher.findAllLeafChildrenWithStructures(chebisearcher.findChebiIdByName("tricarboxylic acid").chebi_id)

ids = [r.chebi_id for r in res]

names = [r.chebi_name for r in res]
print(names)
~~~~

And looking for all molecular entities that have a particular role:
~~~~

res = chebisearcher.findAllWithRole(chebisearcher.findChebiIdByName("nicotinic antagonist").chebi_id)
ids = [r.chebi_id for r in res]

names = [r.chebi_name for r in res]
print(names)

~~~~
Of course, in addition to the pre-built search functions, you can build up your own queries dynamically:

~~~~
res = chebisearcher.findAllByQueryString(chebisearcher.IS_A+chebisearcher.findChebiIdByName("subatomic particle").chebi_id+chebisearcher.AND+chebisearcher.IS_A+chebisearcher.findChebiIdByName("molecular entity").chebi_id)
[r.chebi_name for r in res]

~~~~


If you don't want to search but just want to access the stored database by chebi id (for example), use: 

~~~~ 

from chebidblite import dblite

dbchebi = dblite.ChebiDbLite()
dbchebi.initialize()

chebi_id = "CHEBI:15377"
water = dbchebi.getEntity(chebi_id)
print(water.chebi_name)
# Direct is_a relationships only:
print(water.is_a)
print([dbchebi.getEntity(e).chebi_name for e in water.is_a])

# All recursively populated ancestor IDs for each entity are stored in a separate map
ancestors_of_water = dbchebi.ancestor_map[chebi_id]
print([dbchebi.getEntity(e).chebi_name for e in ancestors_of_water])

~~~~

 


