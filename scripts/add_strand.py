'''
Add strand information to summary file from a genepred source 

Usage:
    python add_strand.py -i <summary_file> -g <genepred_source> -o <output_file>

'''

import argparse
import pandas as pd

def main(args):
    summary = pd.read_csv(args.input, sep='\t')
    summary.columns = ['name', 'chr', 'start', 'end', 'count', 'sum', 'min', 'max', 'mean', 'median', 'std', 'avg_exon_mappabbility', 'link']

    genepred = pd.read_csv(args.genepred, sep='\t', header=None)
    genepred.columns = ['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds']

    summary = summary.merge(genepred[['name', 'strand']], on='name', how='left')
    summary.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add strand information to summary file from a genepred source')
    parser.add_argument('-i', '--input', help='Input summary file', required=True)
    parser.add_argument('-g', '--genepred', help='Genepred source', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)
    args = parser.parse_args()
    main(args)