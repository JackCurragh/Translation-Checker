'''
main script for executing this workflow
'''
import sys
import os

# get path to this scripts directory
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f'{dir_path}/scripts/')

from exclude_unmappable_regions import main as exclude_unmappable_regions
from signal_checker import main as signal_checker
from ranges_csv_to_tx_summary import main as ranges_csv_to_tx_summary
from UCSC_BED_to_ranges_csv import main as UCSC_BED_to_ranges_csv

import argparse
import subprocess
from rich import print

def read_genomic_ranges(genomic_ranges_path: str) -> list:
    '''
    read a genomic ranges csv file

    input:
        genomic_ranges_path: str
            a csv file of genomic ranges

    output:
        genomic_ranges: list
            a list of genomic ranges
    '''
    return pd.read_csv(genomic_ranges_path, sep='\t', header=None).values.tolist()

def split_bed_by_chromosome(bed_file: str, output_dir: str) -> list:
    '''
    create chromosome specific bed files to break up the analysis into smaller chunks
    
    input:
        bed_file: str
            a bed file
        output_dir: str
            a directory to write the output files

    output:
        bed_files: list
            a list of bed files
    '''
    print(f'[bold green]Splitting bed file by chromosome[/bold green]')
    bed_files = []
    with open(bed_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            chromosome = line[0]
            output_file = os.path.join(output_dir, f'{chromosome}.bed')
            with open(output_file, 'a') as o:
                o.write('\t'.join(line) + '\n')
            bed_files.append(output_file)
    print(f'[bold green]Done - {len(bed_files)} files created[/bold green]')
    return bed_files

def convert_bed_to_ranges_csv(bed_file: str, output_file: str) -> str:
    '''
    convert a bed file to a ranges csv file

    input:
        bed_file: str
            a bed file
        output_dir: str
            a directory to write the output files

    output:
        ranges_csv: str
            a ranges csv file
    '''
    print(f'[bold green]Converting bed file to ranges csv file[/bold green]')

    print(f'[bold green]Converting {bed_file} to {output_file}[/bold green]')
    UCSC_BED_to_ranges_csv(bed_file, output_file)
    print(f'[bold green]Done - {output_file} created[/bold green]')
    return output_file


def create_mappability_file(genomic_ranges_path: str, mappability_bw:str, output_dir: str) -> str:
    '''
    create a bigWig file of mappability

    input:
        bed_file: str
            a bed file
        output_dir: str
            a directory to write the output files

    output:
        mappability_file: str
            a bigWig file of mappability
    '''
    print(f'[bold green]Creating mappability file[/bold green]')
    mappability = signal_checker(genomic_ranges_path, mappability_bw, output_dir)
    return mappability


def exclude_unmappable_regions_from_ranges_csv(ranges_csv: str, mappability_file: str, output_dir: str) -> str:
    '''
    exclude unmappable regions from a ranges csv file

    input:
        ranges_csv: str
            a ranges csv file
        mappability_file: str
            a bigWig file of mappability
        output_dir: str
            a directory to write the output files

    output:
        ranges_csv: str
            a ranges csv file
    '''
    print(f'[bold green]Excluding unmappable regions from ranges csv file[/bold green]')
    output_file = os.path.join(output_dir, f'{os.path.basename(ranges_csv).replace(".csv", "")}_no_unmappable_regions.csv')
    exclude_unmappable_regions(mappability_file, output_file)
    print(f'[bold green]Done - {output_file} created[/bold green]')
    return output_file


def create_translation_support_csv(ranges_csv: str, signal_bw: str, output_dir: str) -> str:
    '''
    create a csv file of translation support

    input:
        ranges_csv: str
            a ranges csv file
        signal_bw: str
            a bigWig file of signal
        output_dir: str
            a directory to write the output files

    output:
        translation_support_csv: str
            a csv file of translation support
    '''
    print(f'[bold green]Creating translation support csv file[/bold green]')
    output_file = os.path.join(output_dir, f'{os.path.basename(ranges_csv).replace(".csv", "")}_translation_support.csv')
    output_file = output_dir
    signal_checker(ranges_csv, signal_bw, output_dir, output_file)
    print(output_dir)

    print(f'[bold green]Done - {output_file} created[/bold green]')
    return output_file

def convert_ranges_csv_to_tx_summary(ranges_csv: str, output_dir: str) -> str:
    '''
    convert a ranges csv file to a tx summary file

    input:
        ranges_csv: str
            a ranges csv file
        output_dir: str
            a directory to write the output files

    output:
        tx_summary: str
            a tx summary file
    '''
    print(f'[bold green]Converting ranges csv file to tx summary file[/bold green]')
    output_file = os.path.join(output_dir, f'{os.path.basename(ranges_csv).replace(".csv", "")}_tx_summary.csv')
    ranges_csv_to_tx_summary(ranges_csv, output_file)
    print(f'[bold green]Done - {output_file} created[/bold green]')
    return output_file



def main(args: argparse.Namespace):
    '''
    Parse arguments and execute the workflow
    '''
    # ensure output directory exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # # split the bed file by chromosome
    # bed_files = split_bed_by_chromosome(args.bed, args.output)

    if args.bed:
        print(f'[bold blue]Creating ranges csv file[/bold blue]')
        genomic_ranges = convert_bed_to_ranges_csv(args.bed, f"{args.output}/genomic_ranges.csv")
        print(f'[bold green]Done - {genomic_ranges} created[/bold green]')

    elif args.genomic_ranges:
        print(args.genomic_ranges)
        genomic_ranges = args.genomic_ranges

    # print(f'[bold blue]Creating mappability file[/bold blue]')
    # mappability = create_mappability_file(genomic_ranges, args.mappability, f"{args.output}/mappability.csv")
    # print(f'[bold green]Done - {mappability} created[/bold green]')
    
    # print(f'[bold blue]Excluding unmappable regions from ranges csv file[/bold blue]')
    # ranges_csv = exclude_unmappable_regions_from_ranges_csv(genomic_ranges, mappability, f"{args.output}/mappable_ranges.csv")
    # print(f'[bold green]Done - {ranges_csv} created[/bold green]')
    ranges_csv = genomic_ranges


    print(f'[bold blue]Creating translation support csv file[/bold blue]')
    translation_support = create_translation_support_csv(ranges_csv, args.signal, f"{args.output}/translation_support.csv")
    print(f'[bold green]Done - {translation_support} created[/bold green]')


    print(f'[bold blue]Converting ranges csv file to tx summary file[/bold blue]')
    tx_summary = convert_ranges_csv_to_tx_summary(translation_support, f"{args.output}/tx_summary.csv")
    print(f'[bold green]Done - {tx_summary} created[/bold green]')

    print(f'[bold green]Done - all files created[/bold green]')
    
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Python script to remove unmappable regions from a genomic ranges csv file")
    parser.add_argument("-b", "--bed", type=str, help="A bed annotation file as per UCSC standards")
    parser.add_argument("-g", "--genomic_ranges", type=str, help="A csv file of genomic ranges")
    parser.add_argument("-m", "--mappability", type=str, help="A bigWig file of umap position mappability")
    parser.add_argument("-s", "--signal", type=str, help="A bigWig file of signal")
    parser.add_argument("-o", "--output", type=str, help="A directory path to write output files")
    parser.add_argument("-M", "--mappability_threshold", type=float, help="A threshold for mappability")
    parser.add_argument("-S", "--signal_threshold", type=float, help="A threshold for signal")
    parser.add_argument("-t", "--threads", type=int, help="Number of threads to use")
    args = parser.parse_args()
    main(args)
