import argparse
import os
import copy
from Bio import SeqIO, Phylo
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

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

def load_orthologs(ortho_files):
    orthologs = defaultdict(set)
    ortholog_set = set()
    for ortho_file in ortho_files:
        with open(ortho_file, 'r') as f:
            for line in f:
                cols = line.strip().split('\t')
                orthologs[cols[0]].add(cols[1])
                ortholog_set.add(cols[0])
    return orthologs, ortholog_set

def get_phylip_names(filepath, name_mapping_file):
    basename = os.path.basename(filepath)
    raw_name = os.path.splitext(basename)[0]

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

def process_entry(entry_id, entry_seq, focal_name, orthologs, neighbor_seqs, neighbor_names, tree):
    if entry_id in orthologs:
        output_fasta = f"{entry_id}_n_orthologs.fna"
        output_tree = f"{entry_id}.nwk"
        with open(output_fasta, 'w') as out_f:
            # Write the focal sequence with the FASTA file name in the header
            entry_seq.id = f"{focal_name}"
            entry_seq.description = ""
            entry_seq.seq = entry_seq.seq[:-3] # Remove the terminal stop codon
            linear_seq = str(entry_seq.seq)
            out_f.write(f">{entry_seq.id}\n{linear_seq}\n")

            # Write the ortholog sequences with the FASTA file name in the header
            ortholog_names = {focal_name}
            for ortho_id in orthologs[entry_id]:
                for neighbor_fasta, seqs in neighbor_seqs.items():
                    if ortho_id in seqs:
                        seq = seqs[ortho_id]
                        seq.id = f"{neighbor_names[neighbor_fasta]}"
                        seq.description = ""
                        seq.seq = seq.seq[:-3] # Remove the terminal stop codon
                        linear_seq = str(seq.seq)
                        out_f.write(f">{seq.id}\n{linear_seq}\n")
                        ortholog_names.add(seq.id)

        # Extract the subset of the tree
        subset_tree = extract_subtree(tree, ortholog_names)
        Phylo.write(subset_tree, output_tree, "newick")

def main():
    print("get_alignment_fastas.py\n")
    
    num_cpus = len(os.sched_getaffinity(0))
    print(f"Using {num_cpus} CPUs")

    args = parse_arguments()

    # Compute the basename of the focal FASTA file without the last extension
    focal_name = get_phylip_names(args.fna_a, args.name_mapping)

    # Compute the basenames of the neighbor FASTA files without the last extension
    neighbor_names = {fasta: get_phylip_names(fasta, args.name_mapping) for fasta in args.fna_b}

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
                    neighbor_seqs[fasta] = sequences
            except Exception as exc:
                print(f"{fasta} generated an exception: {exc}")

    print("Sequences loaded")

    # Load orthologs from ortho files
    print("Loading orthologs...")
    orthologs, ortholog_set = load_orthologs(args.ortho)
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
        for entry_id, entry_seq in focal_seqs.items():
            if entry_id in ortholog_set:
                futures.append(executor.submit(process_entry, entry_id, entry_seq, focal_name, orthologs, neighbor_seqs, neighbor_names, tree))

        for future in as_completed(futures):
            future.result()
            processed_count += 1
            if processed_count % 100 == 0 or processed_count == total_count:
                print(f"Processed {processed_count}/{total_count} sequences")

if __name__ == "__main__":
    main()