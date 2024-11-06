#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 16:40:15 2022

@author: paul.roginski
"""


import argparse
import dendropy
import pandas as pd
import sys


# Arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("-names", required=True, help="csv file matching names (col1) and phylip_names (col2)", default=False)
parser.add_argument("-tree", required=True, help="newick file for the input phylogeny tree")
parser.add_argument("-out", required=True, help="output newick file")
args = parser.parse_args()
    
# newick file for the phylogeny tree
tree_file = args.tree

names = pd.read_csv(args.names, names = ["name", "phylip_name"])

tree = dendropy.Tree.get(path=tree_file, schema='newick', preserve_underscores=True)

# Remove internal nodes :
for node in tree :
    node.label=None

for name in names["name"] :
    
    try :
        node_to_change = tree.find_node_with_taxon_label(name)
        node_to_change.taxon.label = names.loc[names['name'] == name, 'phylip_name'].item()
    except AttributeError :
        sys.exit("Taxon " + "'" + name + "'" + " not found in the tree")


# print(tree)
# tree.write(path="/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/SCER_NCBI/lol", schema="newick")
tree.write(path=args.out, schema="newick")