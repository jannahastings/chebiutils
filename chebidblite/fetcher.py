import os
import urllib.request
import gzip

class ChebiFetcher:
    
    CHEBI_OBO_NIGHTLY = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/nightly/chebi.obo.gz'
    CHEBI_OBO_LOCAL = 'chebi.obo'

    def __init__(self):
        self.cacheDir = os.getenv('CHEBIDBLITECACHE', '~')
        self.debug=True
        
    def fetchChebiOBONightly(self):
        if self.debug: print("Fetching ChEBI OBO file from FTP")
        os.makedirs(self.cacheDir,exist_ok=True)
        response = urllib.request.urlopen(self.CHEBI_OBO_NIGHTLY)
        filedata = response.read()
        s_out = gzip.decompress(filedata)
        with open(self.cacheDir+self.CHEBI_OBO_LOCAL, 'wb') as f:
            f.write(s_out)
        if self.debug: print("ChEBI OBO downloaded")




