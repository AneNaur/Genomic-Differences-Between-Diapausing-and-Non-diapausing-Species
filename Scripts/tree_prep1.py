from Bio import AlignIO
from ete3 import Tree
from Bio import Phylo
from Bio import SeqIO
from Bio.AlignIO import PhylipIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from io import StringIO
import re
import os

#NEXUS tree to Newick tree without branch lengths
def nexus_to_newick(nexus_tree, output_tree):
    #load nexus file tree
    tree = Phylo.read(nexus_tree, "nexus")

    #convert nexus to newick for easy processing
    handle = StringIO() 
    Phylo.write(tree, handle, "newick") #creates file like object in memory without saving to disk

    #load/create tree object with ete3
    t = Tree(handle.getvalue())

    #write tree to newick file
    t.write(format=9, outfile=output_tree) #format tells what whould be included in the tree, 9 means only leaf names and no branch lengths


#Paths to data
base = os.getcwd()
in_tree = os.path.join(base, "omm_v12_190tax_CDS_tree.nexus")

#Convert full tree
nexus_to_newick(in_tree, "full_omm_tree.nwk")
