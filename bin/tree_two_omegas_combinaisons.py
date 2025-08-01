#!/usr/bin/env python3

import argparse
import os
from ete3 import Tree

def label_single_node(tree, node_name, mode):
    if mode == "branch":
        suffix =" #1"
    elif mode == "subtree":
        suffix = " $1"
    for node in tree.traverse():
        if node.name == node_name:
            node.name += suffix
    return tree

def main():
    parser = argparse.ArgumentParser(description='Process a Newick tree file.')
    parser.add_argument('file', type=str, help='The Newick tree file')
    parser.add_argument('--ascii', action='store_true', help='Print the ASCII plot')

    args = parser.parse_args()

    # Read the Newick tree
    original_tree = Tree(args.file)
    
    # Get basename without extension for output files
    basename = os.path.splitext(os.path.basename(args.file))[0]

    # Assign unique names to each node
    for i, node in enumerate(original_tree.traverse()):
        if not node.is_leaf():
            node.name = "node" +str(i)

    # Initialize iterator for file naming
    file_counter = 1

    # For each node in the original tree, create a copy of the tree and label the node as "#1"
    for node in original_tree.traverse():
        if not node.is_root():
            for mode in ["branch", "subtree"]:
                if not node.is_leaf() and not node.is_root() or mode != 'subtree':
                    new_tree = original_tree.copy()
                    new_tree = label_single_node(new_tree, node.name, mode)
                    
                    # Write tree to file
                    output_filename = f"{basename}_2omega_{file_counter}.nwk"
                    with open(output_filename, 'w') as f:
                        f.write(new_tree.write(format=1))
                    print(f"Written: {output_filename}")
                    
                    file_counter += 1
                    
                    if args.ascii:
                        print(new_tree.get_ascii(show_internal=True))

if __name__ == "__main__":
    main()