#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import copy
from Bio import SeqIO, Phylo
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys
import logging

logging.basicConfig(level=logging.WARNING)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Build FASTA files with orthologs for each entry in fna_a.")
    parser.add_argument('--fna_a', required=True, help='Path to the focal FASTA file')
    parser.add_argument('--fna_b', required=True, nargs='+', help='Paths to the neighbor FASTA files')
    parser.add_argument('--ortho', required=True, nargs='+', help='Paths to the orthologs TSV files')
    parser.add_argument('--name_mapping', required=True, help='Path to the name mapping TSV file')
    parser.add_argument('--tree', required=True, help='Path to the input Newick tree file')
    return parser.parse_args()

def load_sequences(fasta_file):
    return {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}

def load_orthologs(ortho_files, ortho_files_to_names):
    orthologs = defaultdict(dict)
    ortholog_set = set()
    for ortho_file in ortho_files:
        neighbor = ortho_files_to_names[ortho_file]
        with open(ortho_file, 'r') as f:
            for line in f:
                cols = line.strip().split('\t')
                orthologs[cols[0]][neighbor] = cols[1]
                ortholog_set.add(cols[0])
    return orthologs, ortholog_set

def get_phylip_names(filepath, name_mapping_file, type):

    basename = os.path.basename(filepath)
    raw_name = os.path.splitext(basename)[0]

    if type == 'ortho_file':
        # ex : ortho_pairs/Pfal_iORF_vs_Plasmodium_sp._DRC_iORF_orthologs.tsv
        raw_name = raw_name.split('_vs_')[1]
        raw_name = raw_name.split('_orthologs')[0]

    with open(name_mapping_file, 'r') as f:
        for line in f:
            cols = line.strip().split(',')
            if cols[0] == raw_name:
                return cols[1]

    raise ValueError(f"PHYLIP name for {raw_name} not found in {name_mapping_file}")

def extract_subtree(tree, leaf_names):
    # Create a deep copy of the tree to avoid modifying the original tree
    subtree = copy.deepcopy(tree)

    # Prune the tree to include only the specified leaves
    for leaf in subtree.get_terminals():
        if leaf.name not in leaf_names:
            subtree.prune(leaf)

    return subtree

def process_entry(entry_id, entry_seq, focal_name, orthologs, neighbor_seqs, tree):
    output_fasta = f"{entry_id}_n_orthologs.fna"
    output_tree = f"{entry_id}.nwk"
    
    sequences = {}
    ortholog_names = {focal_name} # for the tree

    # Add the focal sequence
    sequences[focal_name] = str(entry_seq.seq[:-3]) # Remove the terminal stop codon
    
    # Add the ortholog sequences
    for neighbor in orthologs[entry_id]:
        ortho_id = orthologs[entry_id][neighbor]
        try:
            seq = neighbor_seqs[neighbor][ortho_id]
            seq.seq = seq.seq[:-3]
            sequences[neighbor] = str(seq.seq)
            ortholog_names.add(neighbor)
        except KeyError:
            logging.warning(f"Entry ID {ortho_id} not found in the fasta of neighbor {neighbor}.")
            sys.exit(1)        
    
    # Check if all sequences are identical
    if len(set(sequences.values())) == 1:
        logging.warning("All sequences in the output FASTA are identical. The file will not be written.")
        return
    
    # Write the FASTA file
    with open(output_fasta, 'w') as out_f:
        for header, sequence in sequences.items():
            out_f.write(f">{header}\n{sequence}\n")
    
    # Extract the subset of the tree
    subset_tree = extract_subtree(tree, ortholog_names)
    Phylo.write(subset_tree, output_tree, "newick")

def main():
    print("get_alignment_fastas.py\n")
    
    num_cpus = len(os.sched_getaffinity(0))
    print(f"Using {num_cpus} CPUs")

    args = parse_arguments()

    # Compute the basename of the focal FASTA file without the last extension
    focal_name = get_phylip_names(args.fna_a, args.name_mapping, 'fasta')

    # Compute the PHYLIP names for the neighbor FASTA files
    neighbor_fasta_to_names = {fasta: get_phylip_names(fasta, args.name_mapping, 'fasta') for fasta in args.fna_b}

    # Compute the PHYLIP names for the ortholog files
    ortho_files_to_names = {ortho: get_phylip_names(ortho, args.name_mapping, 'ortho_file') for ortho in args.ortho}

    # Filter out non-FASTA files
    fasta_files = [f for f in [args.fna_a] + args.fna_b if f.endswith(('.fna', '.fasta', '.fa'))]

    # Load sequences from focal and neighbor FASTA files using parallel processing
    print("Loading sequences...")

    with ThreadPoolExecutor(max_workers=num_cpus) as executor:
        future_to_fasta = {executor.submit(load_sequences, fasta): fasta for fasta in fasta_files}
        focal_seqs = None
        neighbor_seqs = {}

        for future in as_completed(future_to_fasta):
            fasta = future_to_fasta[future]
            try:
                sequences = future.result()
                if fasta == args.fna_a:
                    focal_seqs = sequences
                else:
                    neighbor = neighbor_fasta_to_names[fasta]
                    neighbor_seqs[neighbor] = sequences
            except Exception as exc:
                print(f"{fasta} generated an exception: {exc}")

    print("Sequences loaded")

    # Load orthologs from ortho files
    print("Loading orthologs...")
    # Remove A_vs_A_orthologs.tsv from ortho_files if it exists
    with open(args.name_mapping, 'r') as f:
        for line in f:
            cols = line.strip().split(',')
            if cols[1] == focal_name:
                focal_raw_name = cols[0]
                break
    A_vs_A_orthologs = f"{focal_raw_name}_vs_{focal_raw_name}_orthologs.tsv"
    print(f"Checking for {A_vs_A_orthologs} in ortholog files...")
    print(f"Ortholog files: {args.ortho}")
    for ortho_file in args.ortho:
        if os.path.basename(ortho_file) == A_vs_A_orthologs:
            print(f"Removing {A_vs_A_orthologs} from ortholog files as it is not needed.")
            args.ortho.remove(ortho_file)
    orthologs, ortholog_set = load_orthologs(args.ortho, ortho_files_to_names)
    print("Orthologs loaded")

    # Load the input tree
    print("Loading tree...")
    tree = Phylo.read(args.tree, "newick")
    print("Tree loaded")

    # Process each entry in fna_a using parallel processing
    print("Processing entries...")

    processed_count = 0
    total_count = len(ortholog_set)

    with ThreadPoolExecutor(max_workers=num_cpus) as executor:
        futures = []
        for entry_id in ortholog_set:
            try:
                entry_seq = focal_seqs[entry_id]
                # Your existing code to process entry_seq
            except KeyError:
                logging.warning(f"Entry ID {entry_id} not found in focal {args.fna_a}.")
                sys.exit(1)
            futures.append(executor.submit(process_entry, entry_id, entry_seq, focal_name, orthologs, neighbor_seqs, tree))

        for future in as_completed(futures):
            future.result()
            processed_count += 1
            if processed_count % 100 == 0 or processed_count == total_count:
                print(f"Processed {processed_count}/{total_count} sequences")

if __name__ == "__main__":
    main()