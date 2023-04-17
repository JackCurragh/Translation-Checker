'''
main script for executing this workflow
'''
import sys
import os

# get path to this scripts directory
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f'{dir_path}/scripts/')

from exclude_unmappable_regions import main as exclude_unmappable_regions
from Signal_Checker import main as signal_checker
from ranges_csv_to_tx_summary import main as ranges_csv_to_tx_summary
from UCSC_BED_to_ranges_csv import main as UCSC_BED_to_ranges_csv
from ranges_csv_to_tx_summary import main as ranges_csv_to_UCSC_BED

import argparse
import subprocess
from rich import print

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

def main(args: argparse.Namespace):
    '''
    Parse arguments and execute the workflow
    '''
    # ensure output directory exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # split the bed file by chromosome
    bed_files = split_bed_by_chromosome(args.bed, args.output)
    
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Python script to remove unmappable regions from a genomic ranges csv file")
    parser.add_argument("-b", "--bed", type=str, help="A bed annotation file as per UCSC standards")
    parser.add_argument("-m", "--mappability", type=str, help="A bigWig file of umap position mappability")
    parser.add_argument("-s", "--signal", type=str, help="A bigWig file of signal")
    parser.add_argument("-o", "--output", type=str, help="A directory path to write output files")
    parser.add_argument("-M", "--mappability_threshold", type=float, help="A threshold for mappability")
    parser.add_argument("-S", "--signal_threshold", type=float, help="A threshold for signal")
    parser.add_argument("-t", "--threads", type=int, help="Number of threads to use")
    args = parser.parse_args()
    main(args)
