#!/usr/bin/env python3

import argparse
import pandas as pd
from scipy.stats import chi2

def main():
    parser = argparse.ArgumentParser(description='Process a codeml stats TSV file.')
    parser.add_argument('file', type=str, help='The codeml stats TSV file')

    args = parser.parse_args()

    # Read the TSV file
    df = pd.read_csv(args.file, sep='\t')

    # Compute the AIC column
    df['AIC'] = 2 * df['np'] - 2 * df['lnL']

    # Initialize the deltalnL and p-value columns
    df['test'] = 'NA'
    df['deltalnL'  ] = 'NA'
    df['p-value'   ] = 'NA'

    # For each model except the 1omega model
    for i in range(1, df.shape[0]):

        # If the model name include "omega_fixed", compare it with the previous model
        if 'omega_fixed' in df.loc[i, 'model']:
            df.loc[i, 'test'] = df.loc[i-1, 'model'] + " better than " + df.loc[i, 'model']
            if df.loc[i, 'lnL'] < df.loc[i-1, 'lnL']:
                df.loc[i, 'deltalnL'] = 2 * (df.loc[i-1, 'lnL'] - df.loc[i, 'lnL'])
                df.loc[i, 'p-value'] = chi2.sf(df.loc[i, 'deltalnL'], df.loc[i-1, 'np'] - df.loc[i, 'np'])
        else:
            # Compare it with the reference 1omega model
            df.loc[i, 'test'] = df.loc[i, 'model'] + " better than " + df.loc[0, 'model']
            if df.loc[i, 'lnL'] > df.loc[0, 'lnL']:
                df.loc[i, 'deltalnL'] = 2 * (df.loc[i, 'lnL'] - df.loc[0, 'lnL'])
                df.loc[i, 'p-value'] = chi2.sf(df.loc[i, 'deltalnL'], df.loc[i, 'np'] - df.loc[0, 'np'])
    
    # Write the result back to the original TSV file
    df.to_csv(args.file, sep='\t', index=False)

if __name__ == "__main__":
    main()