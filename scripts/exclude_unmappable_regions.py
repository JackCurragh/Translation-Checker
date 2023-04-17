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
    region_mappability = pd.read_csv(region_mappability_file, header=1, names=['name','chr', 'start', 'end', 'sum', 'score'])
    return region_mappability

def exclude_unmappable_regions(region_mappability: pd.DataFrame, threshold: float=0.5) -> pd.DataFrame:
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
    mappable_regions = region_mappability[region_mappability['score'] >= threshold]
    print(mappable_regions)
    return mappable_regions

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
    region_mappability = read_region_mappability(args.region_mappability)
    genomic_ranges = exclude_unmappable_regions(region_mappability, threshold=0.5)
    write_genomic_ranges(genomic_ranges, args.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Python script to remove unmappable regions from a genomic ranges csv file')
    parser.add_argument('-m', '--region_mappability', type=str, help='A list of genomic ranges (chr, start, end, mappability)')
    parser.add_argument('-o', '--output', type=str, help='A list of genomic ranges (name, chr, start, end) with unmappable regions removed')
    args = parser.parse_args()
    main(args)