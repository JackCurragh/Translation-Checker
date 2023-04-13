'''
Python script to remove unmappable regions from a genomic ranges csv file

Parameters
----------
genomic_ranges : str
    A list of genomic ranges (name, chr, start, end)
region_mappability : str
    A list of genomic ranges (chr, start, end, sum, score)
output_file : str
    A list of genomic ranges (name, chr, start, end) with unmappable regions removed

Returns
-------
genomic_ranges : pd.DataFrame
    A dataframe of genomic ranges (name chr, start, end)

'''
import argparse
import pandas as pd

def read_genomic_ranges(genomic_ranges_file: str) -> pd.DataFrame:
    '''
    Read a list of genomic ranges (name, chr, start, end)

    Parameters
    ----------
    genomic_ranges_file : str
        A list of genomic ranges (name, chr, start, end)

    Returns
    -------
    genomic_ranges : pd.DataFrame
        A dataframe of genomic ranges (name chr, start, end)
    '''
    genomic_ranges = pd.read_csv(genomic_ranges_file, sep='\t', header=None, names=['name', 'chr', 'start', 'end'])
    return genomic_ranges

def read_region_mappability(region_mappability_file: str) -> pd.DataFrame:
    '''
    Read a list of genomic ranges (chr, start, end, mappability)

    Parameters
    ----------
    region_mappability_file : str
        A list of genomic ranges (chr, start, end, sum, score)

    Returns
    -------
    region_mappability : pd.DataFrame
        A dataframe of genomic ranges (chr, start, end, sum, score)
    '''
    region_mappability = pd.read_csv(region_mappability_file, header=None, names=['chr', 'start', 'end', 'sum', 'score'])
    return region_mappability

def exclude_unmappable_regions(genomic_ranges: pd.DataFrame, region_mappability: pd.DataFrame, threshold: float=0) -> pd.DataFrame:
    '''
    Exclude unmappable regions from a list of genomic ranges

    Parameters
    ----------
    genomic_ranges : pd.DataFrame
        A dataframe of genomic ranges (name chr, start, end)
    region_mappability : pd.DataFrame
        A dataframe of genomic ranges (chr, start, end, mappability)

    Returns
    -------
    genomic_ranges : pd.DataFrame
        A dataframe of genomic ranges (name chr, start, end) with unmappable regions removed
    '''
    genomic_ranges = genomic_ranges.merge(region_mappability, on=['chr', 'start', 'end'], how='left')
    genomic_ranges = genomic_ranges[genomic_ranges['score'] > threshold]
    genomic_ranges = genomic_ranges.drop(columns=['score', 'sum'])
    return genomic_ranges

def write_genomic_ranges(genomic_ranges: pd.DataFrame, output_file: str) -> None:
    '''
    Write a list of genomic ranges (name, chr, start, end) with unmappable regions removed

    Parameters
    ----------
    genomic_ranges : pd.DataFrame
        A dataframe of genomic ranges (name chr, start, end) with unmappable regions removed
    output_file : str
        A list of genomic ranges (name, chr, start, end) with unmappable regions removed
    '''
    genomic_ranges.to_csv(output_file, sep='\t', header=False, index=False)

def main(args):
    genomic_ranges = read_genomic_ranges(args.genomic_ranges)
    print(genomic_ranges)
    region_mappability = read_region_mappability(args.region_mappability)
    print(region_mappability)
    genomic_ranges = exclude_unmappable_regions(genomic_ranges, region_mappability, threshold=0.5)
    print(genomic_ranges)
    write_genomic_ranges(genomic_ranges, args.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Python script to remove unmappable regions from a genomic ranges csv file')
    parser.add_argument('-g', '--genomic_ranges', type=str, help='A list of genomic ranges (name, chr, start, end)')
    parser.add_argument('-m', '--region_mappability', type=str, help='A list of genomic ranges (chr, start, end, mappability)')
    parser.add_argument('-o', '--output', type=str, help='A list of genomic ranges (name, chr, start, end) with unmappable regions removed')
    args = parser.parse_args()
    main(args)