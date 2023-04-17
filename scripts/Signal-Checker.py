"""
Take a list of genomic ranges and return a score for how much signal support there is for each range.

Input:
    - A list of genomic ranges (name, chr, start, end)
    - A Ribo-Seq File (BED or BigWig)

Output:
    - A list of genomic ranges with a score for how much signal support there is for each range. (chr, start, end, score)

Usage:
    python Signal-Checker.py -g <genomic ranges file> -r <ribo-seq file> -o <output file> 

Options:
    -g, --genomic-ranges-file <genomic ranges file>  A list of genomic ranges (chr, start, end)
    -r, --ribo-seq-file <ribo-seq file>             A Ribo-Seq File (BED or BigWig)
    -o, --output-file <output file>                 A list of genomic ranges with a score for how much signal support there is for each range. (chr, start, end, score)
    -h, --help                                      Show this screen.
"""

import argparse
import polars as pl
import pyBigWig
import pandas as pd
import numpy as np

import dask.dataframe as dd
from dask.distributed import Client



def read_genomic_ranges(genomic_ranges_file: str) -> dd.DataFrame:
    """
    Read a list of genomic ranges (name, chr, start, end)

    Parameters
    ----------
    genomic_ranges_file : str
        A list of genomic ranges (name, chr, start, end)

    Returns
    -------
    genomic_ranges : dd.DataFrame
        A dataframe of genomic ranges (name chr, start, end)
    """
    genomic_ranges = dd.read_csv(
        genomic_ranges_file,
        sep="\t",
        header=None,
        names=["name", "chr", "start", "end"],
    )
    return genomic_ranges

def read_ribo_seq_bed(ribo_seq_file: str) -> dd.DataFrame:
    """
    Read a Ribo-Seq File that is in BED format (chr, start, end, score)

    Parameters
    ----------
    ribo_seq_file : str
        A Ribo-Seq File (BED)

    Returns
    -------
    ribo_seq : dd.DataFrame
        A dataframe of ribo-seq data (chr, start, end, score)
    """
    try:
        ribo_seq = dd.read_csv(
            ribo_seq_file,
            sep="\t", 
            header=None, 
            names=["chr", "start", "end", "score"]
        )
        return ribo_seq
        
    except:
        raise Exception("Ribo-Seq file does not appear to be a BED file")


def read_ribo_seq_bigwig(ribo_seq_file: str) -> pyBigWig:
    """
    Read a Ribo-Seq File that is in BigWig format

    Parameters
    ----------
    ribo_seq_file : str
        A Ribo-Seq File (BigWig)

    Returns
    -------
    ribo_seq : pyBigWig
        A pyBigWig object of ribo-seq data
    """
    bw = pyBigWig.open(ribo_seq_file)

    if bw.isBigWig():
        return bw
    else:
        raise Exception("The file is not a BigWig file")


def calculate_signal_support_bed(
    genomic_ranges: pl.DataFrame, ribo_seq: pl.DataFrame
) -> pl.DataFrame:
    """
    Calculate the signal support for each genomic range

    Parameters
    ----------
    genomic_ranges : pl.DataFrame
        A dataframe of genomic ranges (chr, start, end)
    ribo_seq : pl.DataFrame
        A dataframe of ribo-seq data (chr, start, end, score)   

    Returns
    -------
    signal_support : pl.DataFrame
        A dataframe of genomic ranges with a score for how much signal support there is for each range. (chr, start, end, score)
    """
    signal_support = pl.DataFrame(schema={"name": pl.Utf8, "chr": pl.Utf8, "start": pl.Int32, "end": pl.Int32, "sum": pl.Float32, "score": pl.Float32})
    signal_support = {"name": [], "chr": [], "start": [], "end": [], "sum": [], "score": []}
    for row in genomic_ranges.iter_rows(named=True):
        genomic_range = ribo_seq.filter(
            (pl.col('chr') == row["chr"])
            & (pl.col('start')  >= row["start"])
            & (pl.col('end')  <= row["end"])
        )
        if genomic_range.collect().shape[0] == 0:
            signal_support["name"].append(row["name"])
            signal_support["chr"].append(row["chr"])
            signal_support["start"].append(row["start"])
            signal_support["end"].append(row["end"])
            signal_support["sum"].append(0)
            signal_support["score"].append(0)
        else:
            signal_support["name"].append(row["name"])
            signal_support["chr"].append(row["chr"])
            signal_support["start"].append(row["start"])
            signal_support["end"].append(row["end"])
            signal_support["sum"].append(genomic_range.sum().collect()['score'][0])
            signal_support["score"].append(genomic_range.sum().collect()['score'][0] / (row["end"] - row["start"]))

    signal_support = pl.DataFrame(signal_support)
    return signal_support


def calculate_signal_support_bigwig(
    genome_ranges: pd.DataFrame, bw: pyBigWig, cutoff: int = 50
) -> pd.DataFrame:
    """
    Calculate signal support for each genomic range

    Parameters
    ----------
    genomic_ranges : pd.DataFrame
        A dataframe of genomic ranges (chr, start, end)
    bw : pyBigWig
        A pyBigWig object of ribo-seq data

    Returns
    -------
    signal_support : pd.DataFrame
        A dataframe of genomic ranges with a score for how much signal support there is for each range. (chr, start, end, score)
    """
    def generate_entry(row):
        if row['end'] - row['start'] <= cutoff:
            return pd.DataFrame({
                "name": [None],
                "chr": [None],
                "start": [None],
                "end": [None],
                "sum": [None],
                "score": [None],
            })        
        if row['chr'] in bw.chroms():
            try:
                genomic_range = bw.stats(row['chr'], row['start'], row['end'], type="sum")
            except:
                genomic_range = [0]
        else:
            genomic_range = [0]
        if genomic_range[0] is None:
            genomic_range[0] = 0
        return pd.DataFrame({
                "name": [row['name']],
                "chr": [row['chr']],
                "start": [row['start']],
                "end": [row['end']],
                "sum": [genomic_range[0]],
                "score": [genomic_range[0] / (row['end'] - row['start'])],
            })
    signal_support_list = list(genome_ranges.apply(generate_entry, axis=1, meta=[("name", str), ('chr', str), ('start', int), ('end', int), ('sum', float), ('score', float)]))
    signal_support = pd.concat(signal_support_list, ignore_index=True)

    return signal_support


def write_output_pl(signal_support: pl.DataFrame, output_file: str):
    """
    Write output to specified filepath
    """
    signal_support.write_csv(output_file, has_header=False)

def write_output_pd(signal_support: dd.DataFrame, output_file: str):
    """
    Write output to specified filepath
    """
    signal_support.to_csv(output_file, single_file=True, index=False)


def low_memory(args):
    ''' 
    Run checker with low memory usage
    '''
    def check_bigwig_alignment(group, bigwig_path):
        with pyBigWig.open(bigwig_path) as bw:
            results = []
            for _, row in group.iterrows():
                name, chr, start, end, sum, mappability = row
                values = bw.values(chr, int(start), int(end))
                cleaned_values = [x for x in values if str(x) != 'nan']
                sum = np.sum(cleaned_values)
                score = sum / (end - start)
                results.append((name, chr, int(start), int(end), sum, score, mappability))
            return pd.DataFrame(results, columns=["name", "chr", "start", "end", "sum", 'score', 'mappability'])

    # Read the CSV file using Dask
    csv_file = args.genomic_ranges_file
    bigwig_file = args.ribo_seq_file
    ddf = dd.read_csv(csv_file, sep='\t', header=None, names=["name", "chr", "start", "end", "sum", "score"])

    # Group the Dask DataFrame by 'name'
    grouped_ddf = ddf.groupby("name")

    # Iterate through the unique names and process one group at a time
    unique_names = grouped_ddf.name.unique().compute().tolist()

    for idx, name in enumerate(unique_names):
        if idx % 100 == 0:
            print(f"Processing {name} ({idx + 1}/{len(unique_names)})")
        group = grouped_ddf.get_group(name[0])
        result = check_bigwig_alignment(group, bigwig_file)
        with open(args.output_file, "a") as f:
            result.to_csv(f, header=False, index=False)

def main(genomic_ranges_file: str, ribo_seq_file: str, output_file: str, cutoff: int = 50, format: str = None):
    """
    Run checker 
    """
    genomic_ranges = read_genomic_ranges(genomic_ranges_file)

    if not format:
        file_is_bw = pyBigWig.open(ribo_seq_file).isBigWig()
    else:
        file_is_bw = format == "bw"

    if file_is_bw:
        ribo_seq = read_ribo_seq_bigwig(ribo_seq_file)
        signal_support = calculate_signal_support_bigwig(
            genomic_ranges, ribo_seq, cutoff
        )
        write_output_pd(signal_support, output_file)

    else:
        ribo_seq = read_ribo_seq_bed(ribo_seq_file)
        signal_support = calculate_signal_support_bed(
            genomic_ranges, ribo_seq
        )
        write_output_pl(signal_support, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Take a list of genomic ranges and return a score for how much signal support there is for each range."
    )
    parser.add_argument(
        "-g",
        "--genomic-ranges-file",
        type=str,
        required=True,
        help="A list of genomic ranges (chr, start, end)",
    )
    parser.add_argument(
        "-r",
        "--ribo-seq-file",
        type=str,
        required=True,
        help="A Ribo-Seq File (BED or BigWig)",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        type=str,
        required=True,
        help="A list of genomic ranges with a score for how much signal support there is for each range. (chr, start, end, score)",
    )
    parser.add_argument(
        "-f",
        "--format",
        type=str,
        required=False,
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        type=int,
        required=False,
        default=50,
        help="Minimum length of genomic range to calculate signal support",
    )
    parser.add_argument(
        "-l",
        "--low-memory",
        action="store_true",
        help="Use low memory mode",
    )
    args = parser.parse_args()
    if args.low_memory:
        low_memory(args)
    else:
        main(args.genomic_ranges_file, args.ribo_seq_file, args.output_file, args.cutoff, args.format)


    # # Concatenate the results and save to a new CSV file
    # final_df = pd.concat(all_results)
    # final_df.to_csv("cls_annotation/chr21_test.csv", index=False)
