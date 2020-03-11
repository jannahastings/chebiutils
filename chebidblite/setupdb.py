import os
from .fetcher import ChebiFetcher
from .dblite import ChebiDbLite
from .indexer import ChebiIndexer

def prepareCacheAndIndex():
    f = ChebiFetcher()
    f.fetchChebiOBONightly()
    datab = ChebiDbLite()
    datab.initialize(from_cache=False)
    ind = ChebiIndexer()
    ind.buildIndex()
    # all good to go!



# from chebidblite import dblite
# datab = dblite.ChebiDbLite()
# datab.initialize(from_cache=False)