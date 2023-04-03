'''
add extra info from a csv file to the tx summary file.
Join done on tx name 

Parameters
----------
tx_summary_csv : str
    A CSV file of scores for transcripts (name, sum, min, max, mean, median, std)
extra_info_csv : str    
    A CSV file of extra info for transcripts (name, extra_info)

Returns
-------
tx_summary : pd.DataFrame
    A dataframe of scores for transcripts (name, sum, min, max, mean, median, std, extra_info)

    https://gwips.ucc.ie/cgi-bin/hgTracks?db=hg38&position=chr1%3A22564330-22604436
'''

import argparse
import pandas as pd

def read_tx_summary_csv(tx_summary_csv: str) -> pd.DataFrame:
    '''
    Read a CSV file of scores for transcripts (name, sum, min, max, mean, median, std)
    Parameters
    ----------
    tx_summary_csv : str
        A CSV file of scores for transcripts (name, sum, min, max, mean, median, std)
    Returns
    -------
    tx_summary : pd.DataFrame
        A dataframe of scores for transcripts (name, sum, min, max, mean, median, std)
    '''
    tx_summary = pd.read_csv(
        tx_summary_csv,
        sep="\t",
        header=0,
        names=["name", "sum", "min", "max", "mean", "median", "std"],
    )
    return tx_summary

def read_extra_info_csv(extra_info_csv: str) -> pd.DataFrame:
    '''
    Read a CSV file of extra info for transcripts (name, extra_info)
    Parameters
    ----------
    extra_info_csv : str
        A CSV file of extra info for transcripts (name, extra_info)
    Returns
    -------
    extra_info : pd.DataFrame
        A dataframe of extra info for transcripts (name, extra_info)
    '''
    extra_info = pd.read_csv(
        extra_info_csv,
        header=0,
    )
    return extra_info

def generate_gwips_link(tx_summary: pd.DataFrame) -> pd.DataFrame:
    '''
    Generate a link to the GWIPS browser for each transcript
    Parameters
    ----------
    tx_summary : pd.DataFrame
        A dataframe of scores for transcripts (name, sum, min, max, mean, median, std, extra_info)
    Returns
    -------
    tx_summary : pd.DataFrame
        A dataframe of scores for transcripts (name, sum, min, max, mean, median, std, extra_info, gwips_link)
    '''
    tx_summary['gwips_link'] = tx_summary['name'].apply(lambda x: f'https://gwips.ucc.ie/cgi-bin/hgTracks?db=hg38&position={x}')
    return tx_summary


def add_extra_info(tx_summary: pd.DataFrame, extra_info: pd.DataFrame) -> pd.DataFrame:
    '''
    Add extra info to tx summary file
    Join done on tx name
    Parameters
    ----------
    tx_summary : pd.DataFrame
        A dataframe of scores for transcripts (name, sum, min, max, mean, median, std)
    extra_info : pd.DataFrame
        A dataframe of extra info for transcripts (name, extra_info)
    Returns
    -------
    tx_summary : pd.DataFrame
        A dataframe of scores for transcripts (name, sum, min, max, mean, median, std, extra_info)
    '''
    new_columns = ['exon_number']
    for col in new_columns:
        tx_summary[col] = None

    for i, row in extra_info.iterrows():

        tx_summary.loc[tx_summary['name'] == row[3], 'exon_number'] = row[4]
    return tx_summary

def write_tx_summary_csv(tx_summary: pd.DataFrame, tx_summary_csv: str) -> None:
    '''
    Write a CSV file of scores for transcripts (name, sum, min, max, mean, median, std, extra_info)
    Parameters
    ----------
    tx_summary : pd.DataFrame
        A dataframe of scores for transcripts (name, sum, min, max, mean, median, std, extra_info)
    tx_summary_csv : str
        A CSV file of scores for transcripts (name, sum, min, max, mean, median, std, extra_info)
    '''
    tx_summary.to_csv(tx_summary_csv, index=False)

def main():
    parser = argparse.ArgumentParser(description='add extra info from a csv file to the tx summary file. Join done on tx name')
    parser.add_argument('--tx_summary_csv', type=str, help='A CSV file of scores for transcripts (name, sum, min, max, mean, median, std)')
    parser.add_argument('--extra_info_csv', type=str, help='A CSV file of extra info for transcripts (name, extra_info)')
    parser.add_argument('--tx_summary_csv_out', type=str, help='A CSV file of scores for transcripts (name, sum, min, max, mean, median, std, extra_info)')
    args = parser.parse_args()

    tx_summary = read_tx_summary_csv(args.tx_summary_csv)
    extra_info = read_extra_info_csv(args.extra_info_csv)
    tx_summary = add_extra_info(tx_summary, extra_info)
    write_tx_summary_csv(tx_summary, args.tx_summary_csv_out)

if __name__ == "__main__":
    main()