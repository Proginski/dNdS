#!/usr/bin/env python3

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Process a codeml stats TSV file.')
    parser.add_argument('file', type=str, help='The codeml stats TSV file')
    parser.add_argument('--pvalue', type=float, default=0.05, help='The p-value threshold')

    args = parser.parse_args()

    # Read the TSV file
    df = pd.read_csv(args.file, sep='\t', keep_default_na=False)

    # Add a new "selection" column initialized with "NA"
    df['selection'] = 'NA'

    for index, row in df.iterrows():
        if row['p-value'] != "NA" and float(row['p-value']) <= args.pvalue:
            
            # For models with "1omega_fixed1" or "free_omega" in their name
            if '1omega_fixed1' in row['model'] or 'free_omega' in row['model']:
                    df.loc[index, 'selection'] = 'true'

            # For models with "_2omega" and no "_fixed" in their name
            elif '_2omega' in row['model'] and '_fixed' not in row['model']:
                fixed_model_name = row['model'].rsplit('_', 1)[0] + '_fixed1_' + row['model'].rsplit('_', 1)[1]
                fixed_model = df[df['model'] == fixed_model_name]
                # Test the p-value of the fixed model
                if fixed_model['p-value'].values[0] != "NA" and float(fixed_model['p-value'].values[0]) <= args.pvalue:
                    # Set selection to "true"
                    df.at[index, 'selection'] = 'true'

    # Write the result back to the original TSV file
    df.to_csv(args.file, sep='\t', index=False)

if __name__ == "__main__":
    main()