#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO
import time

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate FASTA or TSV files for orthologs.")
    parser.add_argument('--fna_a', required=True, help='Path to the focal FASTA file')
    parser.add_argument('--fna_b', required=True, nargs='+', help='Paths to the neighbor FASTA files')
    parser.add_argument('--ortho', required=True, nargs='+', help='Paths to the orthologs TSV files')
    # parser.add_argument('--output_type', required=True, choices=['fasta', 'tsv'], help='Output type: fasta or tsv')
    # parser.add_argument('--output_dir', required=True, help='Directory to save the output files')
    parser.add_argument('--output', required=True, help='Output TSV file')
    return parser.parse_args()

def load_sequences(fasta_file):
    return {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}

def load_orthologs(ortho_files):
    orthologs = {}
    for ortho_file in ortho_files:
        with open(ortho_file, 'r') as f:
            for line in f:
                focal_orf, neighbor_orf = line.strip().split('\t')
                if focal_orf not in orthologs:
                    orthologs[focal_orf] = []
                orthologs[focal_orf].append((ortho_file, neighbor_orf))
    return orthologs

# def generate_fasta_output(focal_seqs, neighbor_seqs, orthologs, output_dir):
#     for focal_orf, neighbors in orthologs.items():
#         output_fasta = os.path.join(output_dir, f"{focal_orf}.fasta")
#         with open(output_fasta, 'w') as out_f:
#             SeqIO.write(focal_seqs[focal_orf], out_f, "fasta")
#             for ortho_file, neighbor_orf in neighbors:
#                 neighbor_fasta = neighbor_seqs[ortho_file]
#                 SeqIO.write(neighbor_fasta[neighbor_orf], out_f, "fasta")

def generate_tsv_output(focal_seqs, neighbor_seqs, orthologs, output_file):
    with open(output_file, 'w') as out_f:
        for focal_orf, neighbors in orthologs.items():
            line = [focal_orf, str(focal_seqs[focal_orf].seq)]
            for ortho_file, neighbor_orf in neighbors:
                neighbor_fasta = neighbor_seqs[ortho_file]
                line.extend([neighbor_orf, str(neighbor_fasta[neighbor_orf].seq)])
            out_f.write('\t'.join(line) + '\n')

def main():
    global_start = time.time()
    print("get_alignment_fastas.py")

    args = parse_arguments()

    # Load sequences
    start = time.time()
    print("Loading sequences...", end=' ')
    focal_seqs = load_sequences(args.fna_a)
    neighbor_seqs = {os.path.basename(fasta): load_sequences(fasta) for fasta in args.fna_b}
    print(f"Done ({time.time() - start:.2f}s)")

    # Load orthologs
    start = time.time()
    print("Loading orthologs...", end=' ')
    orthologs = load_orthologs(args.ortho)
    print(f"Done ({time.time() - start:.2f}s)")

    # Generate output
    start = time.time()
    print("Generating output...", end=' ')
    # if args.output_type == 'fasta':
    #     os.makedirs(args.output_dir, exist_ok=True)
    #     generate_fasta_output(focal_seqs, neighbor_seqs, orthologs, args.output_dir)
    # elif args.output_type == 'tsv':
    #     output_file = os.path.join(args.output_dir, "output.tsv")
    #     generate_tsv_output(focal_seqs, neighbor_seqs, orthologs, output_file)
    generate_tsv_output(focal_seqs, neighbor_seqs, orthologs, args.output)
    print(f"Done ({time.time() - start:.2f}s)")

    print(f"Total time: {time.time() - global_start:.2f}s")

if __name__ == "__main__":
    main()