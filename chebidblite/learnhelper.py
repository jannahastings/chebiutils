import os
import pandas as pd
from chebidblite import searcher
from chebidblite import dblite
from rdkit import Chem
import pickle
import networkx as nx
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from colour import Color
import numpy
import random
from ontoutils import robot_wrapper
from pronto import Ontology

# Clean chebi.obo for Pronto (in case of syntax errors):
#sed '/PDBeChem/d' chebi.obo > chebi-tmp.obo
#sed '/MetaCyc/d' chebi-tmp.obo > chebi.obo



# ChebiDataPreparer prepares ChEBI data for machine learning
# FIRST TIME:
# dprep = ChebiDataPreparer(load_from_cache=False)
# dprep.buildDefaultDataMaps()
# AFTERWARDS:
# dprep = ChebiDataPreparer()
#
# chemdata = dprep.getDataForLearning(class_size=100,number_to_select=100)
#
class ChebiDataPreparer:

    def __init__(self,load_from_cache=True,debug=False):
        self.chebisearcher = searcher.ChebiSearcher()
        self.db = dblite.ChebiDbLite()

        # Indexes and maps
        self.classes_to_members={}  # Map of parent classes to their leaf structure holder classes
        self.members_to_classes={}  # Map of leaf structure chebis to their parent classes
        self.chebis_to_smiles={}    # Map of leaf structures to smiles classes
        self.bad_chebis_to_smiles = {}  # Map of chebi Ids and smiles for which fingerprints cannot be generated
        self.smiles_to_fps={}       # Map of smiles structures to fingerprints
        self.debug=debug

        if load_from_cache:
            self.__loadFromCache__()

    def __saveToCache__(self):
        if self.debug: print("Saving pickled data for ChebiDataPreparer")
        if os.path.exists(self.db.cacheDir+"/LearningDataPreparer/"):
            shutil.rmtree((self.db.cacheDir+"/LearningDataPreparer/"))
        os.makedirs(self.db.cacheDir+"/LearningDataPreparer/")
        with open(self.db.cacheDir+"/LearningDataPreparer/classes_to_members.pkl",'wb') as output:
            pickle.dump(self.classes_to_members,output)
        with open(self.db.cacheDir+"/LearningDataPreparer/members_to_classes.pkl",'wb') as output:
            pickle.dump(self.members_to_classes,output)
        with open(self.db.cacheDir+"/LearningDataPreparer/chebis_to_smiles.pkl",'wb') as output:
            pickle.dump(self.chebis_to_smiles,output)
        with open(self.db.cacheDir+"/LearningDataPreparer/bad_chebis_to_smiles.pkl",'wb') as output:
            pickle.dump(self.bad_chebis_to_smiles,output)
        with open(self.db.cacheDir+"/LearningDataPreparer/smiles_to_fps.pkl",'wb') as output:
            pickle.dump(self.smiles_to_fps,output)


    def __loadFromCache__(self):
        if self.debug: print("Loading pickled data for ChebiDataPreparer")
        with open(self.db.cacheDir+"/LearningDataPreparer/classes_to_members.pkl",'rb') as input:
            self.classes_to_members = pickle.load(input)
        with open(self.db.cacheDir+"/LearningDataPreparer/members_to_classes.pkl",'rb') as input:
            self.members_to_classes = pickle.load(input)
        with open(self.db.cacheDir+"/LearningDataPreparer/chebis_to_smiles.pkl",'rb') as input:
            self.chebis_to_smiles = pickle.load(input)
        with open(self.db.cacheDir+"/LearningDataPreparer/bad_chebis_to_smiles.pkl",'rb') as input:
            self.bad_chebis_to_smiles = pickle.load(input)
        with open(self.db.cacheDir+"/LearningDataPreparer/smiles_to_fps.pkl",'rb') as input:
            self.smiles_to_fps = pickle.load(input)


    def buildDefaultDataMaps(self,rootname="molecular entity",fpSize=1024):
        self.prepareClassToStructureMap(rootname)
        self.generateStructureFingerprints(fpSize)
        self.__saveToCache__()

    def prepareClassToStructureMap(self, rootname="molecular entity"):
        # From ChEBI 'dblite' get all the non-leaf rootname classes that have leaf nodes with structures. Get the structures as SMILES.
        molents = self.chebisearcher.findAllByQueryString(self.chebisearcher.IS_A+ \
                    self.chebisearcher.findChebiIdByName(rootname).chebi_id+ \
                    self.chebisearcher.AND+"leaf_node:False")
        self.classes_to_members[rootname] = set()
        for m in molents:
            str_ch = self.chebisearcher.findAllLeafChildrenWithStructures(m.chebi_id)
            if str_ch is not None:
                self.classes_to_members[m.chebi_id]=set([s.chebi_id for s in str_ch])
                if self.debug:
                    print(m.chebi_id,"has no. of structures:",len(str_ch))
                for s in str_ch:
                    self.chebis_to_smiles[s.chebi_id] = s.smiles
                    if s.chebi_id not in self.members_to_classes:
                        self.members_to_classes[s.chebi_id] = set()
                    self.members_to_classes[s.chebi_id].add(m.chebi_id)
                    self.classes_to_members[rootname].add(m.chebi_id)

    # Generate the fingerprints for the whole once only for each chebiId from the smiles
    # Warning: slow
    def generateStructureFingerprints(self,fpSize=1024):
        for chem in self.chebis_to_smiles:
            smiles = self.chebis_to_smiles[chem]
            if smiles not in self.smiles_to_fps:
                ms = Chem.MolFromSmiles(smiles)
                if ms is not None:
                    fp = Chem.RDKFingerprint(ms, fpSize=fpSize)
                    self.smiles_to_fps[smiles]= list(fp.ToBitString())
                else:
                    self.bad_chebis_to_smiles[chem] = smiles



    # TODO: THE SELECTION NEEDS REWORKING
    # Need functionality to rapidly select N randomly selected classes for testing
    # structure-based learning, at a lower or higher level, with X members, maximally
    # disjoint, no shared structures between the classes.
    #
    # Something like a distance metric between classes based on the count of shared
    # leaves (structures) would be very helpful.
    #
    # Then build the data frame by under-sampling favouring the most specific classes
    # Order from smallest to largest
    # Keep a note of the ones not used aside for testing and validation.
    def getDataForLearning(self,class_size,number_to_select):

        # Work from the smallest classes to the biggest classes
        clslens = {k:len(v) for k,v in self.classes_to_members.items() if len(v)>class_size}
        clsorder = sorted(clslens.items() ,  key=lambda x: x[1])

        used_children = set()
        skipped_classes = {}

        chemdata = pd.DataFrame(columns=numpy.arange(1027))
        rowindex = 0
        clscount = 0
        for elem in clsorder:
            parent_id = elem[0]
            children = [c for c in self.classes_to_members[parent_id] if c not in self.bad_chebis_to_smiles]
            children_to_use = [c for c in children if c not in used_children]
            if len(children_to_use) < class_size:
                #n_missing = class_size - len(children_to_use)
                print("Not enough never-seen children in class ",parent_id, ", only ",len(children_to_use),", SKIPPING")
                skipped_classes[parent_id]=len(children_to_use)
                #if n_missing > len(children)-len(children_to_use):
                #    n_missing = len(children)-len(children_to_use)
                #duplicate_children = random.sample([c for c in children if c not in children_to_use],n_missing)
                #children_to_use = children_to_use + duplicate_children
                # SKIP THIS CLASS, DON'T ADD DUPLICATES
                continue
            else:
                if len(children_to_use) > class_size:
                    # Randomly select class_size of them
                    children_to_use = random.sample(children_to_use,class_size)
                # Count up the classes
                clscount = clscount + 1
                # Add them to the data frame
                for c in children_to_use:
                    smiles = self.chebis_to_smiles[c]
                    fp = self.smiles_to_fps[smiles]
                    chemdata.loc[rowindex] = [parent_id, c, smiles] + fp
                    rowindex = rowindex + 1
                    if rowindex % 1000 == 0:
                        print("Processing row",rowindex, "---",parent_id)
                    # Now remember that we used these ones
                    used_children.add(c)
                if clscount == number_to_select:
                    return( chemdata )  # Favour more specific classes. Stop when found enough
        return ( chemdata )


    # Gets a table of rows (molecules) x columns (classes) in which every
    # class has at least 'class_size' members,  and every molecule
    # is flagged in all the classes it belongs to.
    #
    # The number of classes parameter is the initial seed, but
    # additional classes may need
    # to be added in order to complete a spanning tree of all the molecules
    def getDataForDeepLearning(self, class_size, number_to_select):
        # Work from the smallest classes to the biggest classes, as the smaller
        # classes are more specific thus more likely to be learnable
        clslens = {k:len(v) for k,v in self.classes_to_members.items() if len(v)>class_size}
        clsorder = sorted(clslens.items() ,  key=lambda x: x[1])
        selected_classes = [c for c,_v in clsorder[:number_to_select] ]

        # Get these classes together with their
        # (randomly selected from the total) class_size members.
        all_members = set().union(*[random.sample(self.classes_to_members[c],class_size) for c in selected_classes])


        chemdata = pd.DataFrame(columns=['MOLECULEID','SMILES']+selected_classes)

        rowindex = 0

        for m in all_members:
            smiles = self.chebis_to_smiles[m]
            in_classes = [True if m in self.classes_to_members[c] else False for c in selected_classes]
            chemdata.loc[rowindex] = [m, smiles] + in_classes
            rowindex = rowindex + 1
            if rowindex % 1000 == 0:
                print("Processing row",rowindex, "---",m)

        return chemdata



class ChebiOntologySubsetter:
# Extract a subset of ChEBI that is for just the
# classes that are included in the experiment, for visualising
# the F1 score as a colour on a network plot.
    def __init__(self):
        self.chebislim = None
        self.classes_in = None

    def createSubsetFor(self,classes_in):
        self.classes_in = set(classes_in)
        print(len(classes_in))

        with open("classes_in.txt", 'w') as outfile:
            for c in self.classes_in:
                outfile.writelines(c+"\n")

        rw = robot_wrapper.RobotWrapper(robotcmd='/Users/hastingj/Work/Onto/robot/robot')

#get_ontology_cmd = 'curl -L http://purl.obolibrary.org/obo/chebi.obo > chebi.obo'
#rw.__executeCommand__(get_ontology_cmd)

        extract_cmd = [rw.robotcmd,  "extract --method MIREOT ",
                "--input chebi.obo", "--lower-terms classes_in.txt",
                "--intermediates minimal","--output chebi-slim.obo"]

        rw.__executeCommand__(" ".join(extract_cmd))


        self.chebislim = Ontology("chebi-slim.obo")

    # For the ontology subset, print a network image for those classes
    # including colours for some numeric value assigned to classes
    # colour_nums = dictionary of chebi IDs to some number between 0 and 1
    # num_cols = how many colours to generate in the range from colour start to end
    def printImageOfSubset(self,image_name, colour_nums, num_cols=100,colour_start="red",colour_end="green"):
        red = Color(colour_start)
        colors = list(red.range_to(Color(colour_end),num_cols))

        G=nx.DiGraph()

        for term in self.chebislim.terms():
            chebi_id = term.id
            chebi_name = term.name
            parents = set([t.id for rel in term.relationships for t in term.relationships[rel]  if rel.name == 'is a'])
            definition = term.definition
            color = 'grey'
            if chebi_id in self.classes_in:
                cnum = colour_nums[chebi_id] #f1val = classifres[chebi_id]['f1-score']
                color = colors[int(num_cols*cnum)-1]
            G.add_node(chebi_name,color=color)
            for p in parents:
                G.add_edge(chebi_name,self.chebislim[p].name)

        pdot = nx.drawing.nx_pydot.to_pydot(G)

        for i, node in enumerate(pdot.get_nodes()):
            node.set_shape('box')
            node.set_fontcolor('black')
            #node.set_fillcolor(node.color)
            node.set_style('rounded, filled')
            #node.set_color(node.color)

        png_path = image_name #"chebi-slim-vis.png"
        pdot.write_png(png_path,prog='/usr/local/bin/dot')


