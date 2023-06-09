'''
Script to convert a bed file outputted from GWIPS-viz or UCSC genome browser to a csv file of genomic ranges

Usage:
    python3 bed_to_ranges_csv.py -i <bed file> -o <output file>

Options:
    -i  --bed-file <bed file>           A bed file
    -o, --output-file <output file>      A bed file
    -h, --help                           Show this screen.

'''

import argparse
import pandas as pd
import numpy as np

def read_bed(bed: str) -> pd.DataFrame:
    '''
    Read a bed file

    Parameters
    ----------
    bed : str
        A bed file

    Returns
    -------
    bed : pd.DataFrame
        A dataframe of bed data
    '''
    column_names = ['chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
    bed = pd.read_csv(bed, sep='\t', header=None, names=column_names)
    return bed

def bed_to_csv(bed: pd.DataFrame) -> pd.DataFrame:
    '''
    Convert a bed dataframe to a csv dataframe of genomic ranges

    Parameters
    ----------
    bed : pd.DataFrame
        A dataframe of bed data

    Returns
    -------
    bed : pd.DataFrame
        A dataframe of bed data for genomic coding regions
    '''
    def generate_ranges(row):
        start = row['start'] + np.array([x for x in str(row['blockStarts']).split(',') if x], dtype=int)
        end = start + np.array([x for x in str(row['blockSizes']).split(',') if x], dtype=int)
        return pd.DataFrame({
            'name': [row['name']] * len(start),
            'chr': [row['chr']] * len(start),
            'start': start,
            'end': end
        })

    genomic_ranges_list = bed.apply(generate_ranges, axis=1).tolist()
    genomic_ranges = pd.concat(genomic_ranges_list, ignore_index=True)

    return genomic_ranges


def write_bed(bed: pd.DataFrame, output_file: str):
    '''
    Write a bed dataframe to a file

    Parameters
    ----------
    bed : pd.DataFrame
        A dataframe of bed data
    output_file : str
        A bed file
    '''
    bed.to_csv(output_file, sep='\t', header=False, index=False)

def main(input_file: str, output_file: str):
    bed = read_bed(input_file)
    bed = bed_to_csv(bed)
    write_bed(bed, output_file)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to convert a genePred file to a bed file')
    parser.add_argument('-i', '--input', help='A genePred file', required=True)
    parser.add_argument('-o', '--output', help='A bed file', required=True)
    args = parser.parse_args()
    main(args.input, args.output)

