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
import polars as pl
import pyBigWig


def read_genomic_ranges(genomic_ranges_file: str) -> pl.DataFrame:
    """
    Read a list of genomic ranges (name, chr, start, end)

    Parameters
    ----------
    genomic_ranges_file : str
        A list of genomic ranges (name, chr, start, end)

    Returns
    -------
    genomic_ranges : pl.DataFrame
        A dataframe of genomic ranges (name chr, start, end)
    """
    genomic_ranges = pl.read_csv(
        genomic_ranges_file,
        sep="\t",
        has_header=False,
        new_columns=["name", "chr", "start", "end"],
    )
    return genomic_ranges


def read_ribo_seq_bed(ribo_seq_file: str) -> pl.DataFrame:
    """
    Read a Ribo-Seq File that is in BED format (chr, start, end, score)

    Parameters
    ----------
    ribo_seq_file : str
        A Ribo-Seq File (BED)

    Returns
    -------
    ribo_seq : pl.DataFrame
        A dataframe of ribo-seq data (chr, start, end, score)
    """
    try:
        ribo_seq = pl.scan_csv(
            ribo_seq_file,
            sep="\t", 
            has_header=False, 
            dtypes={"chr": pl.Utf8, "start": pl.Int32, "end": pl.Int32, "score": pl.Float32}
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
    genomic_ranges: pl.DataFrame, ribo_seq: pl.DataFrame
) -> pl.DataFrame:
    """
    Calculate the translation support for each genomic range

    Parameters
    ----------
    genomic_ranges : pl.DataFrame
        A dataframe of genomic ranges (chr, start, end)
    ribo_seq : pl.DataFrame
        A dataframe of ribo-seq data (chr, start, end, score)   

    Returns
    -------
    translation_support : pl.DataFrame
        A dataframe of genomic ranges with a score for how much translation support there is for each range. (chr, start, end, score)
    """
    translation_support = pl.DataFrame(schema={"name": pl.Utf8, "chr": pl.Utf8, "start": pl.Int32, "end": pl.Int32, "sum": pl.Float32, "score": pl.Float32})
    translation_support = {"name": [], "chr": [], "start": [], "end": [], "sum": [], "score": []}
    for row in genomic_ranges.iter_rows(named=True):
        genomic_range = ribo_seq.filter(
            (pl.col('chr') == row["chr"])
            & (pl.col('start')  >= row["start"])
            & (pl.col('end')  <= row["end"])
        )
        if genomic_range.collect().shape[0] == 0:
            translation_support["name"].append(row["name"])
            translation_support["chr"].append(row["chr"])
            translation_support["start"].append(row["start"])
            translation_support["end"].append(row["end"])
            translation_support["sum"].append(0)
            translation_support["score"].append(0)
        else:
            translation_support["name"].append(row["name"])
            translation_support["chr"].append(row["chr"])
            translation_support["start"].append(row["start"])
            translation_support["end"].append(row["end"])
            translation_support["sum"].append(genomic_range.sum().collect()['score'][0])
            translation_support["score"].append(genomic_range.sum().collect()['score'][0] / (row["end"] - row["start"]))

    translation_support = pl.DataFrame(translation_support)
    return translation_support


def calculate_translation_support_bigwig(
    genome_ranges: pl.DataFrame, bw: pyBigWig
) -> pl.DataFrame:
    """
    Calculate translation support for each genomic range

    Parameters
    ----------
    genomic_ranges : pl.DataFrame
        A dataframe of genomic ranges (chr, start, end)
    bw : pyBigWig
        A pyBigWig object of ribo-seq data

    Returns
    -------
    translation_support : pl.DataFrame
        A dataframe of genomic ranges with a score for how much translation support there is for each range. (chr, start, end, score)
    """
    translation_support = pl.DataFrame(columns=["name", "chr", "start", "end", "sum"])
    for index, row in genome_ranges.iterrows():
        genomic_range = bw.stats(row["chr"], row["start"], row["end"], type="sum")
        if genomic_range[0] is None:
            genomic_range[0] = 0
        entry = pl.DataFrame({
                "name": [row["name"]],
                "chr": [row["chr"]],
                "start": [row["start"]],
                "end": [row["end"]],
                "sum": [genomic_range[0]],
                "score": [genomic_range[0] / (row["end"] - row["start"])],
            })
        translation_support = pl.concat(
            [entry, translation_support], axis=0
        )

    return translation_support


def write_output(translation_support: pl.DataFrame, output_file: str):
    """
    Write output to specified filepath
    """
    translation_support.write_csv(output_file, has_header=False)


def main(args: argparse.Namespace):
    """
    Run translation checker 
    """
    genomic_ranges = read_genomic_ranges(args.genomic_ranges_file)

    if not args.format:
        file_is_bw = pyBigWig.open(args.ribo_seq_file).isBigWig()
    else:
        file_is_bw = args.format == "bw"

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
    parser.add_argument(
        "-f",
        "--format",
        type=str,
        required=False,
    )
    args = parser.parse_args()
    main(args)
