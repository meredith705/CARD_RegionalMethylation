import pandas as pd
import numpy as np
from datetime import datetime
import argparse
import sys

"""
    Takes the 2 phased cohort beds after intersection with bed region and combines counts across haplotypes
    
    usage: sumhaplotypes.py [-h] -a INPUT_H1_COHORT_BED -b INPUT_H2_COHORT_BED -o OUTPUT_FILE [-m MIN_COV]

    Combine haploype read coverage and filter positions from the group BED that don't pass coverage or cpg minimum filters.

    optional arguments:
      -h, --help            show this help message and exit
      -a INPUT_H1_COHORT_BED, --input_h1_cohort_bed INPUT_H1_COHORT_BED
                            Path to the hap1 cohort combined input TSV file with header.
      -b INPUT_H2_COHORT_BED, --input_h2_cohort_bed INPUT_H2_COHORT_BED
                            Path to the hap1 cohort combined input TSV file with header.
      -o OUTPUT_FILE, --output_file OUTPUT_FILE
                            output filename.
      -m MIN_COV, --min_cov MIN_COV
                            Minimum coverage; default=5x.

"""


def log_time(message):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}")



def combine_haplotypes(f1, f2, outfile, mincov):    

    log_time(f'load files')
    file1 = pd.read_csv(f1, sep="\t") #, compression="gzip")
    file1 = file1.drop_duplicates()
    file2 = pd.read_csv(f2, sep="\t") #, compression="gzip")
    file2 = file2.drop_duplicates()

    log_time(f'merge files')
    # Merge on chrom, start, end
    bed_cols = ["#chrom", "start", "end"]
    merged = pd.merge(file1, file2, on=bed_cols, suffixes=('_f1', '_f2'))
    print(merged.head())

    log_time(f'make new df')
    # copy to result with first 3 columns
    summed = merged[bed_cols].copy()


    # Select all sample coverage columns from each file
    sample_cols_f1 = merged.filter(regex='GRCh38_1').columns
    sample_cols_f2 = merged.filter(regex='GRCh38_2').columns

    print('columns 1:',sample_cols_f1[:3] )
    print('columns 2:', sample_cols_f2[:3])

    # # Sum each pair of corresponding columns
    # for c1, c2 in zip(sample_cols_f1, sample_cols_f2):
    #     print(c1, c1.split("_"))
    #     new_col = c1.replace('GRCh38_1', '')  # remove suffix for clean column names
    #     summed[new_col] = merged[c1] + merged[c2]

    # Ensure matched pairs of columns for each haplotype file 
    assert len(sample_cols_f1) == len(sample_cols_f2), "Mismatch in sample column counts"

    log_time(f'sum coverages')
    # Sum sample columns
    summed_sample_df = pd.concat(
        [
            (merged[c1] + merged[c2]).rename(c1.replace("GRCh38_1_", ""))
            for c1, c2 in zip(sample_cols_f1, sample_cols_f2)
        ],
        axis=1
    )

    log_time(f'make summed df')
    # Final output DataFrame: coordinate + summed sample data
    summed = pd.concat([merged[bed_cols], summed_sample_df], axis=1)

    log_time(f'calculate modFraction and mask on mincov of {mincov}x')
    # Calculate the modFraction column and mask samples with < min cov
    cols = list(summed.columns)

    for i, col in enumerate(cols):
        if col.endswith("_modFraction"):
            # check correct collumns are present
            if i > 0 and i < len(cols) - 1:
                prev_col = cols[i - 1]
                next_col = cols[i + 1]
                if prev_col.endswith("_validCov") and next_col.endswith("_modReads"):
                    if (prev_col.split("_")[1] == next_col.split("_")[1]) :
                        # Calculate modFraction
                        # summed[col] = round( ((summed[next_col] / summed[prev_col]) * 100), 2 )

                        # Check for mincov filter then calculate modFraction
                        valid_cov = summed[prev_col]
                        mod_reads = summed[next_col]

                        # Mask low coverage regions 
                        masked_mod_values = np.where( valid_cov >= mincov, 
                                                    ( round((mod_reads/valid_cov)*100, 2)),
                                                    "." )

                        summed[col] = masked_mod_values
                        
                        
                    else:
                        print('samples dont match?:',prev_col, next_col, col)


    log_time(f'write out df')
    # Write output to  file
    summed.to_csv(outfile, sep="\t", index=False, header=True, chunksize=100000) #, compression="gzip")


if __name__ == "__main__":

    # Make an argument parser
    parser = argparse.ArgumentParser(description="Combine haploype read coverage and filter positions from the group BED that don't pass coverage or cpg minimum filters.")
    parser.add_argument(
        "-a","--input_h1_cohort_bed",
        type=str,
        required=True,
        help="Path to the hap1 cohort combined input TSV file with header."
    )

    parser.add_argument(
        "-b","--input_h2_cohort_bed",
        type=str,
        required=True,
        help="Path to the hap1 cohort combined input TSV file with header."
    )

    parser.add_argument(
        "-o","--output_file",
        type=str,
        required=True,
        help="output filename."
    )

    parser.add_argument(
        "-m","--min_cov",
        type=int,
        required=False,
        default=5,
        help="Minimum coverage; default=5x."
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Parse arguments
    args = parser.parse_args()

    combine_haplotypes(args.input_h1_cohort_bed, args.input_h2_cohort_bed, args.output_file, int(args.min_cov) )






