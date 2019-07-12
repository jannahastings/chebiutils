
import urllib.request
import gzip

class ChebiFetcher:
    
    CHEBI_OBO_NIGHTLY = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/nightly/chebi.obo.gz'
    CHEBI_OBO_LOCAL = 'chebi.obo'

    def __init__(self):
        pass
        
    def fetchChebiOBONightly(self):
        response = urllib.request.urlopen(self.CHEBI_OBO_NIGHTLY)
        filedata = response.read()
        s_out = gzip.decompress(filedata)
        with open(self.CHEBI_OBO_LOCAL, 'wb') as f:
            f.write(s_out)




