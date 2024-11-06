#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script turns an orthologues files from OrthoFinder into a 2 columns TSV file.

Designed for the dNdS Nextflow pipeline

The input file is tab separated, has a header line and three columns.
The script remove the header and the first column.
Now orthologues can be one to one or one to many or many to many.
So in columns two and three the value can be multiple sequences separated by a comma.
The script consider each sequence from the second column and associate it (write) to ***ONLY THE FIRST*** sequence of the third column.
"""

import argparse
import csv

def process_orthologues(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        next(reader, None)  # skip the headers

        for row in reader:
            sequences_2 = row[1].split(', ')
            sequences_3 = row[2].split(', ')

            for seq_2 in sequences_2:
                writer.writerow([seq_2, sequences_3[0]])

parser = argparse.ArgumentParser(description='Process OrthoFinder orthologues file')
parser.add_argument('input_file', type=str, help='Input OrthoFinder orthologues file')
parser.add_argument('output_file', type=str, help='Output 2-column TSV file')

args = parser.parse_args()

process_orthologues(args.input_file, args.output_file)