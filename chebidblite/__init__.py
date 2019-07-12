__all__ = ["dblite","fetcher","indexer","searcher"]

name = "chebidblite"

import os

if 'CHEBIDBLITECACHE' not in os.environ.keys():
    os.environ['CHEBIDBLITECACHE'] = os.path.expanduser('~/Library/Caches/chebidblite/')
    
    
