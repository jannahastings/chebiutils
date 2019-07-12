import itertools
import os
import shutil
import pickle
from pronto import *
from .fetcher import ChebiFetcher

class ChebiEntity:
    
    def __init__(self, chebi_id, chebi_name = None, definition = None,
                 synonyms = None, 
                 alt_ids = None, xrefs = None,
                 stars = None, has_struc = False, smiles = None, 
                 has_role = None, is_a = None, children = None, has_part = None,
                 is_conjugate_base_of = None, 
                 is_conjugate_acid_of = None, 
                 is_enantiomer_of = None, 
                 has_functional_parent = None,
                 has_parent_hydride = None, 
                 is_substituent_group_from = None):
        self.chebi_id = chebi_id
        self.chebi_name = chebi_name
        self.definition = definition
        self.synonyms = synonyms
        self.alt_ids = alt_ids
        self.xrefs = xrefs
        self.stars = stars
        # entities with structures only. Todo add charge, mass etc.
        self.has_struc = has_struc
        self.smiles = smiles
        # ontology relationships (direct, outward direction only)
        self.has_role = has_role
        self.is_a = is_a
        self.children = children # only the is-a children
        self.has_part = has_part
        self.is_conjugate_base_of = is_conjugate_base_of
        self.is_conjugate_acid_of = is_conjugate_acid_of
        self.is_enantiomer_of = is_enantiomer_of
        self.has_functional_parent = has_functional_parent
        self.has_parent_hydride = has_parent_hydride
        self.is_substituent_group_from = is_substituent_group_from


class ChebiDbLite:
    """This class creates a simple in-memory cache of ChEBI's data"""
    
    PICKLE_DIR = "chebidblitecache"
    
    def __init__(self):
        self.cacheDir = os.getenv('CHEBIDBLITECACHE', '~')
        self.data_dict = {}
        self.ancestor_map = {}
        self.role_map = {}
        self.secondary_id_map = {}
        self.initialized = False
        self.debug = True
        
    """ Read the OBO file and map the data into the cache """
    def _loadChebiFromOBO(self):
        if self.debug: print("Loading ChEBI data from OBO file...")
        onto = Ontology(self.cacheDir+ChebiFetcher.CHEBI_OBO_LOCAL)  #todo: check exists, maybe try download using Fetcher if not, or else sensible message...
        if self.debug: print("Parsing ChEBI data into in-memory cache...") 
        for t in onto.terms:
            term = onto[t]
            chebi_id = term.id
            chebi_name = term.name
            definition = str(term.desc)
            if 'subset' in term.other.keys():
                stars = term.other['subset'][0]
            else:
                stars = None
            if 'alt_id' in term.other.keys():
                alt_ids = set(term.other['alt_id'])
            else:
                alt_ids = None
            if len(term.synonyms)>0:
                synonyms = set()
                for s in term.synonyms:
                    synonyms.add(s.desc)
            else:
                synonyms = None
            if 'xref' in term.other.keys():
                xrefs = set(term.other['xref'])
            else:
                xrefs = None
            # entities with structures only
            if 'property_value' in term.other.keys():
                smiles_str = next((s for s in term.other['property_value'] if 'smiles' in s), None)
            else: 
                smiles_str = None
            has_struc = (smiles_str is not None)
            #print(term,'has_struc',has_struc,'smiles_str',smiles_str)
            if has_struc:
                smiles = smiles_str.split('"')[1::2]
                if len(smiles) ==0:
                    smiles = smiles_str.split(' ')[1::2]
                if len(smiles) > 0:
                    smiles = smiles[0]
            else:
                smiles = None
            # relationships
            relnames = str(term.relations.keys())
            if 'has_role' in relnames:
                roles = [term.relations[rel].id for rel in term.relations if rel.obo_name == 'has_role']
                has_role = set(itertools.chain.from_iterable(roles)) 
            else:
                has_role = None
            if 'has_part' in relnames:
                parts = [term.relations[rel].id for rel in term.relations if rel.obo_name == 'has_part']
                has_part = set(itertools.chain.from_iterable(parts)) 
            else: 
                has_part = None
            if 'is_conjugate_base_of' in relnames:
                conjbases = [term.relations[rel].id for rel in term.relations if rel.obo_name == 'is_conjugate_base_of']
                is_conjugate_base_of = set(itertools.chain.from_iterable(conjbases)) 
            else:
                is_conjugate_base_of = None
            if 'is_conjugate_acid_of' in relnames:
                conjacids = [term.relations[rel].id for rel in term.relations if rel.obo_name == 'is_conjugate_acid_of']
                is_conjugate_acid_of = set(itertools.chain.from_iterable(conjacids)) 
            else:
                is_conjugate_acid_of = None
            if 'is_enantiomer_of' in relnames:
                enants = [term.relations[rel].id for rel in term.relations if rel.obo_name == 'is_enantiomer_of']
                is_enantiomer_of = set(itertools.chain.from_iterable(enants))
            else:
                is_enantiomer_of = None
            if 'has_functional_parent' in relnames:
                funcps = [term.relations[rel].id for rel in term.relations if rel.obo_name == 'has_functional_parent']
                has_functional_parent = set(itertools.chain.from_iterable(funcps))
            else: 
                has_functional_parent = None
            if 'has_parent_hydride' in relnames:
                parhds = [term.relations[rel].id for rel in term.relations if rel.obo_name == 'has_parent_hydride']
                has_parent_hydride = set(itertools.chain.from_iterable(parhds))
            else:
                has_parent_hydride = None
            if 'is_substituent_group_from' in relnames:
                subsgs = [term.relations[rel].id for rel in term.relations if rel.obo_name == 'is_substituent_group_from']
                is_substituent_group_from = set(itertools.chain.from_iterable(subsgs))
            else:
                is_substituent_group_from = None
            # direct is_a parents:
            if 'is_a' in relnames:
                isas = [term.relations[rel].id for rel in term.relations if rel.obo_name == 'is_a']
                is_a = set(itertools.chain.from_iterable(isas))
            else:
                is_a = None
            # direct is_a children:
            if len(term.children)>0:
                children = set(term.children.id)
            else:
                children = None
            # build the entity for the cache
            entity = ChebiEntity(chebi_id = chebi_id, 
                                 chebi_name = chebi_name, 
                                 synonyms = synonyms, 
                                 alt_ids = alt_ids, xrefs = xrefs,
                                 stars = stars, has_struc = has_struc, 
                                 smiles = smiles, has_role = has_role, 
                                 is_a = is_a, children = children, 
                                 has_part = has_part,
                                 is_conjugate_base_of = is_conjugate_base_of,
                                 is_conjugate_acid_of = is_conjugate_acid_of,
                                 is_enantiomer_of = is_enantiomer_of, 
                                 has_functional_parent = has_functional_parent,
                                 has_parent_hydride = has_parent_hydride, 
                                 is_substituent_group_from = is_substituent_group_from)
            # store in self.data_dict with index the main ID
            self.data_dict[t] = entity

    # Recursively expand the ancestor map with parents of parents
    def _expandAncestorMap(self,entry,visited):
        parents = self.ancestor_map[entry]
        if len(parents)>0:
            for p in parents:
                if p not in visited:
                    visited = self._expandAncestorMap(p,visited)
            for p in parents:
                self.ancestor_map[entry] = self.ancestor_map[entry].union(self.ancestor_map[p])
        visited.add(entry)
        return visited
        
    def _buildAncestorMap(self):
        if self.debug: print("Build ancestor map...")     
        # Prepare the ancestor map
        self.ancestor_map = {}
        for t in self.data_dict.keys():
            self.ancestor_map[t] = set()
        # Add direct parents to the ancestor map. 
        for t in self.data_dict.keys():
            entity = self.data_dict[t]
            if entity.is_a:
                for p in entity.is_a:
                    self.ancestor_map[t].add(p)
        visited = set()
        for t in self.data_dict.keys():
            if self.debug:
                if len(visited)%1000 == 0:
                    print("Building ancestor map for entry number: ",len(visited))
            if t not in visited: 
                visited = self._expandAncestorMap(t,visited)

    def _buildRoleMap(self):        
    # Todo test if ancestor map has already been built
        if self.debug: print("Build role map...")
        # Prepare the initial role map (direct roles)
        self.role_map={}
        for t in self.data_dict.keys():
            entity = self.data_dict[t]
            if entity.has_role:
                self.role_map[t] = entity.has_role
            else:
                self.role_map[t] = set()
                
        # Add inherited roles from term's is-a parents
        for t in self.data_dict.keys():
            ancestors = self.ancestor_map[t]
            for a in ancestors:
                if len(self.role_map[a])>0:
                    self.role_map[t] = self.role_map[t].union(self.role_map[a])
        
        # Add inherited roles from role parents for asserted roles
        all_roles = set(itertools.chain.from_iterable(self.role_map.values()))
        for r in all_roles:
            r_ancestors = self.ancestor_map[r]
            where_used = {k:v for (k,v) in self.role_map.items() if r in v}
            for k in where_used.keys():
                self.role_map[k]=self.role_map[k].union(r_ancestors)
        
    """ Secondary id map maps from alt_id back to primary ID"""
    def _buildSecondaryIdMap(self):
        if self.debug: print("Build secondary ID map...")
        for t in self.data_dict.keys():
            entity = self.data_dict[t]
            if entity.alt_ids:
                for alt_id in entity.alt_ids:
                    self.secondary_id_map[alt_id] = t


    def _saveToCache(self):
        if self.debug: print("Saving pickled cache data")
        if os.path.exists(self.cacheDir+self.PICKLE_DIR):
            shutil.rmtree(self.cacheDir+self.PICKLE_DIR)
        os.makedirs(self.cacheDir+self.PICKLE_DIR)
        with open(self.cacheDir+self.PICKLE_DIR+"/data_dict.pkl",'wb') as output:
            pickle.dump(self.data_dict,output)
        with open(self.cacheDir+self.PICKLE_DIR+"/ancestor_map.pkl",'wb') as output:
            pickle.dump(self.ancestor_map,output)
        with open(self.cacheDir+self.PICKLE_DIR+"/role_map.pkl",'wb') as output:    
            pickle.dump(self.role_map,output)
        with open(self.cacheDir+self.PICKLE_DIR+"/secondary_id_map.pkl",'wb') as output:
            pickle.dump(self.secondary_id_map,output)

    def _loadFromCache(self):
        if self.debug: print("Loading pickled cache data")
        with open(self.cacheDir+self.PICKLE_DIR+"/data_dict.pkl",'rb') as input:
            self.data_dict = pickle.load(input)
        with open(self.cacheDir+self.PICKLE_DIR+"/ancestor_map.pkl",'rb') as input:
            self.ancestor_map = pickle.load(input)
        with open(self.cacheDir+self.PICKLE_DIR+"/role_map.pkl",'rb') as input:
            self.role_map = pickle.load(input)
        with open(self.cacheDir+self.PICKLE_DIR+"/secondary_id_map.pkl",'rb') as input:
            self.secondary_id_map = pickle.load(input)
        

    # TODO manage state so that the dblite cannot be used unless it has been initialized
    def initialize(self,from_cache=True):
        if from_cache and os.path.exists(self.cacheDir+self.PICKLE_DIR): # the files are found
            if self.debug: print("Initializing from cache")
            self._loadFromCache()
        else:
            if self.debug: print("Initializing by reloading")
            self._loadChebiFromOBO()
            self._buildAncestorMap()
            self._buildRoleMap()
            self._buildSecondaryIdMap()
            self._saveToCache()
        self.initialized = True
        if self.debug: print ("Initialization completed.")
        
    
        
    # TODO add proper error handling etc. what if not found etc. or not initialized...
    def getEntity(self,chebiId):
        if not self.initialized:
            print("Error! Cannot use ChebiDbLite without initializing first.")
            return None
        if chebiId in self.data_dict.keys():
            return (self.data_dict[chebiId])
        if chebiId in self.secondary_id_map.keys():
            return (self.data_dict[self.secondary_id_map[chebiId]])
        return None
        
