import argparse
import pandas as pd


def overlap(row_a: pd.Series, row_b: pd.Series) -> bool:
    """
    Checks if two genomic ranges overlap.

    Args:
        row_a: A pandas Series representing a row from dataframe A.
        row_b: A pandas Series representing a row from dataframe B.

    Returns:
        A boolean value indicating whether the genomic ranges overlap.
    """
    return (((row_a.end > row_b.start) and (row_a.end < row_b.end)) or ((row_a.start > row_b.start) and (row_a.start < row_b.end)))


def filter_translation_file(translation_file_path: str, genomic_ranges_path: str, filtered_file_path: str, chunk_size: int = 1000) -> None:
    """
    Filters entries from file A whose range overlaps with any feature in file B, and writes the resulting
    dataframe to a new csv file.

    Args:
        translation_file_path: A string representing the path to the summary file for the translation.
        genomic_ranges_path: A string representing the path to file B.
        filtered_file_path: A string representing the path to the output file that will contain the filtered entries.
        chunk_size: An integer representing the number of rows to process at a time.

    Returns:
        None.
    """
    # Read the files in chunks
    translation_file_columns = ["name","chr", "start", "end", "count", "sum", "min", "max", "mean", "median", "std", "avg_exon_mappabbility", "link"]
    genomic_ranges_columns = ['chr', 'start', 'end']

    for idx, chunk_a in enumerate(pd.read_csv(translation_file_path, chunksize=chunk_size, header=None, names=translation_file_columns)):
        # Create an empty list to store the rows that do not overlap with genomic_ranges
        filtered_rows = []
        for idx_a, row_a in chunk_a.iterrows():
            overlap_flag = False
            for chunk_b in pd.read_csv(genomic_ranges_path, chunksize=chunk_size, header=None, names=genomic_ranges_columns):
                for idx_b, row_b in chunk_b.iterrows():
                    if overlap(row_a, row_b):
                        overlap_flag = True
                        break
                if overlap_flag:
                    break
            if not overlap_flag:
                filtered_rows.append(row_a)

        if filtered_rows:
            print(f"Chunk {idx} of size {chunk_size}\t {len(filtered_rows)} rows written")
            # Concatenate the filtered rows from this chunk to the output dataframe
            filtered_chunk = pd.concat(filtered_rows, axis=1)
            # Write the filtered chunk to a new file
            filtered_chunk.to_csv(filtered_file_path, mode='a', header=False, index=False)
        else:
            print(f'Chunk {idx} of size {chunk_size} No rows in chunk overlap with genomic ranges.')


if __name__ == '__main__':
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Filter file A based on overlap with features in file B')
    parser.add_argument('--translation_file', type=str, required=True, help='Path to file A')
    parser.add_argument('--genomic_ranges', type=str, required=True, help='Path to file B')
    parser.add_argument('--filtered_file', type=str, required=True, help='Path to output filtered file')
    parser.add_argument('--chunk_size', type=int, default=1000, help='Number of rows to process at a time')

    # Parse command line arguments
    args = parser.parse_args()

    # Filter file A based on overlap with features in file B
    filter_translation_file(args.translation_file, args.genomic_ranges, args.filtered_file, args.chunk_size)