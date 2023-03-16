"""
Take a list of genomic ranges and return a score for how much translation support there is for each range.

Input:
    - A list of genomic ranges (name, chr, start, end)
    - A Ribo-Seq File (BED or BigWig)

Output:
    - A list of genomic ranges with a score for how much translation support there is for each range. (chr, start, end, score)

Usage:
    python Translation-Checker.py -g <genomic ranges file> -r <ribo-seq file> -o <output file> 

Options:
    -g, --genomic-ranges-file <genomic ranges file>  A list of genomic ranges (chr, start, end)
    -r, --ribo-seq-file <ribo-seq file>             A Ribo-Seq File (BED or BigWig)
    -o, --output-file <output file>                 A list of genomic ranges with a score for how much translation support there is for each range. (chr, start, end, score)
    -h, --help                                      Show this screen.
"""

import argparse
import pandas as pd
import pyBigWig


def read_genomic_ranges(genomic_ranges_file: str) -> pd.DataFrame:
    """
    Read a list of genomic ranges (name, chr, start, end)

    Parameters
    ----------
    genomic_ranges_file : str
        A list of genomic ranges (name, chr, start, end)

    Returns
    -------
    genomic_ranges : pd.DataFrame
        A dataframe of genomic ranges (name chr, start, end)
    """
    genomic_ranges = pd.read_csv(
        genomic_ranges_file,
        sep="\t",
        header=None,
        names=["name", "chr", "start", "end"],
    )
    return genomic_ranges


def read_ribo_seq_bed(ribo_seq_file: str) -> pd.DataFrame:
    """
    Read a Ribo-Seq File that is in BED format (chr, start, end, score)

    Parameters
    ----------
    ribo_seq_file : str
        A Ribo-Seq File (BED)

    Returns
    -------
    ribo_seq : pd.DataFrame
        A dataframe of ribo-seq data (chr, start, end, score)
    """
    try:
        ribo_seq = pd.read_csv(
            ribo_seq_file, sep="\t", header=None, names=["chr", "start", "end", "score"]
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


def calculate_translation_support_bed(
    genomic_ranges: pd.DataFrame, ribo_seq: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate the translation support for each genomic range

    Parameters
    ----------
    genomic_ranges : pd.DataFrame
        A dataframe of genomic ranges (chr, start, end)
    ribo_seq : pd.DataFrame
        A dataframe of ribo-seq data (chr, start, end, score)   

    Returns
    -------
    translation_support : pd.DataFrame
        A dataframe of genomic ranges with a score for how much translation support there is for each range. (chr, start, end, score)
    """
    translation_support = pd.DataFrame(columns=["name", "chr", "start", "end", "score"])
    for index, row in genomic_ranges.iterrows():
        genomic_range = ribo_seq[
            (ribo_seq["chr"] == row["chr"])
            & (ribo_seq["start"] >= row["start"])
            & (ribo_seq["end"] <= row["end"])
        ]
        translation_support = translation_support.append(
            {
                "name": row["name"],
                "chr": row["chr"],
                "start": row["start"],
                "end": row["end"],
                "score": genomic_range["score"].sum(),
            },
            ignore_index=True,
        )
    return translation_support


def calculate_translation_support_bigwig(
    genome_ranges: pd.DataFrame, bw: pyBigWig
) -> pd.DataFrame:
    """
    Calculate translation support for each genomic range

    Parameters
    ----------
    genomic_ranges : pd.DataFrame
        A dataframe of genomic ranges (chr, start, end)
    bw : pyBigWig
        A pyBigWig object of ribo-seq data

    Returns
    -------
    translation_support : pd.DataFrame
        A dataframe of genomic ranges with a score for how much translation support there is for each range. (chr, start, end, score)
    """
    translation_support = pd.DataFrame(columns=["name", "chr", "start", "end", "score"])
    for index, row in genome_ranges.iterrows():
        genomic_range = bw.stats(row["chr"], row["start"], row["end"], type="sum")
        translation_support = translation_support.append(
            {
                "name": row["name"],
                "chr": row["chr"],
                "start": row["start"],
                "end": row["end"],
                "score": genomic_range[0],
            },
            ignore_index=True,
        )
    return translation_support


def write_output(translation_support: pd.DataFrame, output_file: str):
    """
    Write output to specified filepath
    """
    translation_support.to_csv(output_file, sep="\t", header=None, index=False)


def main(args: argparse.Namespace):
    """
    Run translation checker 
    """
    genomic_ranges = read_genomic_ranges(args.genomic_ranges_file)

    file_is_bw = pyBigWig.open(args.ribo_seq_file).isBigWig()

    if file_is_bw:
        ribo_seq = read_ribo_seq_bigwig(args.ribo_seq_file)
        translation_support = calculate_translation_support_bigwig(
            genomic_ranges, ribo_seq
        )
    else:
        ribo_seq = read_ribo_seq_bed(args.ribo_seq_file)
        translation_support = calculate_translation_support_bed(
            genomic_ranges, ribo_seq
        )
    write_output(translation_support, args.output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Take a list of genomic ranges and return a score for how much translation support there is for each range."
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
        help="A list of genomic ranges with a score for how much translation support there is for each range. (chr, start, end, score)",
    )
    args = parser.parse_args()
    main(args)
