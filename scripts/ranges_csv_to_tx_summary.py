'''
Take a CSV file of scores for ranges and output a summary of the scores for each unique transcript

Input:
    - A CSV file of scores for ranges (name, chr, start, end, sum, score)

Output:
    - A CSV file of scores for transcripts (name, sum, min, max, mean, median, std)

Usage:
    python ranges_csv_to_tx_summary.py -i <input file> -o <output file>

'''

import argparse
import pandas as pd

def read_ranges_csv(ranges_csv: str) -> pd.DataFrame:
    '''
    Read a CSV file of scores for ranges (name, chr, start, end, sum, score)

    Parameters
    ----------
    ranges_csv : str
        A CSV file of scores for ranges (name, chr, start, end, sum, score)

    Returns
    -------
    ranges : pd.DataFrame
        A dataframe of scores for ranges (name, chr, start, end, sum, score)
    '''
    ranges = pd.read_csv(
        ranges_csv,
        # sep="\t",
        header=0,
        names=["name", "chr", "start", "end", "sum", "score", "mappability"],
    )
    return ranges

def ranges_csv_to_tx_summary(ranges: pd.DataFrame) -> pd.DataFrame:
    '''
    Take a dataframe of scores for ranges and output a summary of the scores for each unique transcript

    Parameters
    ----------
    ranges : pd.DataFrame
        A dataframe of scores for ranges (name, chr, start, end, sum, score)

    Returns
    -------
    tx_summary : pd.DataFrame
        A dataframe of scores for transcripts (name, chr, min_start, max_stop, sum, min, max, mean, median, std)
    '''
    # Summarize the scores for each transcript
    tx_summary = ranges.groupby("name").agg(
        {
            "chr": "first",
            "start": "min",
            "end": "max",
            "name": "size",
            "sum": "sum",
            "score": ["min", "max", "mean", "median", "std"],
            "mappability": "mean",
        }
    )
    tx_summary.columns = ["chr", "start", "end", "count", "sum", "min", "max", "mean", "median", "std", "avg_exon_mappabbility"]

    return tx_summary



def write_tx_summary_csv(tx_summary: pd.DataFrame, tx_summary_csv: str) -> None:
    '''
    Write a CSV file of scores for transcripts (name, sum, min, max, mean, median, std)

    Parameters
    ----------
    tx_summary : pd.DataFrame
        A dataframe of scores for transcripts (name, sum, min, max, mean, median, std)
    tx_summary_csv : str
        A CSV file of scores for transcripts (name, sum, min, max, mean, median, std)

    Returns
    -------
    None
    '''
    tx_summary.to_csv(tx_summary_csv, sep="\t", header=True, index=True)

def make_gwips_link(tx_summary: pd.DataFrame) -> pd.DataFrame:
    '''
    Make a link to the GWIPS website for each transcript

    Parameters
    ----------
    tx_summary : pd.DataFrame
        A dataframe of scores for transcripts (name, sum, min, max, mean, median, std)

    Returns
    -------
    tx_summary : pd.DataFrame
        A dataframe of scores for transcripts (name, sum, min, max, mean, median, std, link)
    '''
    base_url = "https://gwips.ucc.ie/cgi-bin/hgTracks?db=hg38" #&position=chr6%3A52497408-52577060
    links = []
    for idx, row in tx_summary.iterrows():
        link = f"{base_url}&position={row.chr}%3A{row.start}-{row.end}"
        links.append(link)

    tx_summary["link"] = links
    return tx_summary

def main(input, output):
    """
    Convert a CSV file of scores for ranges to a CSV file of scores for transcripts
    """
    ranges = read_ranges_csv(input)
    tx_summary = ranges_csv_to_tx_summary(ranges)
    # sort by sum
    tx_summary = tx_summary.sort_values(by=["sum"], ascending=True)
    # add a link to the GWIPS website
    tx_summary = make_gwips_link(tx_summary)
    write_tx_summary_csv(tx_summary, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="A CSV file of scores for ranges (name, chr, start, end, sum, score)", required=True)
    parser.add_argument("-o", "--output", help="A CSV file of scores for transcripts (name, sum, min, max, mean, median, std)", required=True)
    args = parser.parse_args()

    main(args.input, args.output)