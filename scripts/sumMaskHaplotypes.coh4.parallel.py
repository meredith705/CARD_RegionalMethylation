#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import sys
import gzip
from datetime import datetime
import argparse
from multiprocessing import Pool
import subprocess
from io import StringIO


"""
    Takes the 2 phased cohort beds after intersection with bed region and combines counts across haplotypes
    
    usage: sumhaplotypes.py [-h] -a INPUT_H1_COHORT_BED -b INPUT_H2_COHORT_BED -c INPUT_UNGROUPED_COHORT_BED -o OUTPUT_FILE [-m MIN_COV]

    Combine haploype read coverage and filter positions from the group BED that don't pass coverage or cpg minimum filters.

    optional arguments:
      -h, --help            show this help message and exit
      -a INPUT_H1_COHORT_BED, --input_h1_cohort_bed INPUT_H1_COHORT_BED
                            Path to the hap1 cohort combined input TSV file with header.
      -b INPUT_H2_COHORT_BED, --input_h2_cohort_bed INPUT_H2_COHORT_BED
                            Path to the hap2 cohort combined input TSV file with header.
      -c INPUT_UNGROUPED_COHORT_BED, --input_ungrouped_cohort_bed INPUT_UNGROUPED_COHORT_BED
                            Path to the ungrouped cohort combined input TSV file with header.
      -o, --output_dir OUTPUT_DIR
                        output directory.
      -n, --cohort_name COHORT_NAME
                            output cohort name for naming output files.
      -l, --chromosomes [CHROMOSOMES ...]
                            List of chromosomes to process (optional). If not provided, all chromosomes will be used.
      -t, --threads THREADS
                            threads; default 3.
      -m, --min_cov MIN_COV
                            Minimum coverage; default=5x.
      -w, --write_croms WRITE_CROMS
                            Boolean to determine if filtered positions per chromosome are written out; default False

"""


def log_time_top(message):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}")

def get_logger(chrom):
    log_file = open(f"{chrom}.log", "a")
    def log_time(message):
        log_file.write(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}" + "\n")
        log_file.flush()
    return log_time, log_file



def combine_haplotypes(chrom, f1, f2, f3, output_dir, cohort, mincov, write_per_chrom):    
    # set up log file per chrom
    log_time, log_file = get_logger(chrom)

    log_time(f'process {chrom}')
    files = [f1, f2, f3]
    input_dataframes = []

    for f in files:
        log_time(f'load file: {f}')

        # query the chromosome from the bgzipped BED files
        cmd = ["tabix", f, chrom]

        # isolate the chrom subset
        try:
            # log_time(f'command {cmd}')
            in_f_chrom = subprocess.check_output(cmd, text=True)
            # log_time(f'zcat {in_tsv}')
            # header = subprocess.check_output(f"zcat {in_tsv} | head -1", shell=True, text=True).strip().split("\t")     ## fix before leaving local 
            header = subprocess.check_output(f"gzip -dc {f} | head -1", shell=True, text=True).strip().split("\t")
        except subprocess.CalledProcessError as e:
            print(f"Failed tsv on {f} {chrom}: {e}")
            return 

        # load the file into a dataframe and store it in the list of dfs
        try:
            log_time(f'load in {f} {chrom} subset')
            fdf = pd.read_csv(StringIO(in_f_chrom), sep='\t')
            fdf.columns = header
            fdf = fdf.drop_duplicates()
            input_dataframes.append(fdf)
        except Exception as e:
            log_time(f"Error reading {f} {chrom} files {e}" )
            return p

    log_time(f'store files')
    file1, file2, file3 = input_dataframes

    # file1 = pd.read_csv(f1, sep="\t") #, compression="gzip")
    # file1 = file1.drop_duplicates()
    # file2 = pd.read_csv(f2, sep="\t") 
    # file2 = file2.drop_duplicates()
    # file3 = pd.read_csv(f3, sep="\t") 
    # file3 = file3.drop_duplicates()

    log_time(f'merge files')
    # Merge on chrom, start, end
    bed_cols = ["#chrom", "start", "end"]
    merged = file1.merge( file2, on=bed_cols, suffixes=('_f1', '_f2')) \
               .merge( file3, on=bed_cols, suffixes=('', '_f3'), how="outer").fillna(0)
    print(merged.head())

    log_time(f'make new df')
    # copy to result with first 3 columns
    summed = merged[bed_cols].copy()


    # Select all sample coverage columns from each file
    sample_cols_f1 = merged.filter(regex='ref_1_').columns
    sample_cols_f2 = merged.filter(regex='ref_2_').columns
    sample_cols_f3 = merged.filter(regex='ref_ungrouped_').columns

    print('columns 1:',sample_cols_f1[:3] )
    print('columns 2:', sample_cols_f2[:3])
    print('columns 3:', sample_cols_f3[:3])


    # Ensure matched pairs of columns for each haplotype file 
    assert len(sample_cols_f1) == len(sample_cols_f2) == len(sample_cols_f3), "Mismatch in sample column counts"

    log_time(f'sum coverages')
    # Sum sample columns
    summed_sample_df = pd.concat(
        [
            (merged[c1] + merged[c2] + merged[c3]).rename(c1.replace("ref_1_", ""))
            for c1, c2, c3 in zip(sample_cols_f1, sample_cols_f2, sample_cols_f3)
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

    if write_per_chrom:
        # write out the chrom dataframe
        outfile = f"{cohort}_mergedMeth_{chrom}.tsv"
        log_time(f'write out df to {output_dir}/{outfile}')
        # Write output to  file
        summed.to_csv(f'{output_dir}/{outfile}', sep="\t", index=False, header=True, chunksize=100000) #, compression="gzip")
    else:
        log_time(f"finished {chrom}")

    return summed

 # process each chromosome in parallel
def process_chromosomes_in_parallel(output_dir, input_h1_cohort_bed, input_h2_cohort_bed, input_ungrouped_cohort_bed, cohort, min_cov, write_croms, chromosomes, nprocesses):

    # create a pool of worker processes
    with Pool(processes=nprocesses) as pool:

        results = pool.starmap(combine_haplotypes, [(chrom, input_h1_cohort_bed, input_h2_cohort_bed, input_ungrouped_cohort_bed, output_dir,cohort, min_cov, write_croms) for chrom in chromosomes])

    # combine results
    log_time(f'process parallel results')

    combined_meth = pd.concat([result for result in results], ignore_index=True)

    log_time(f'write out combined data')
    # write out the filtered dataframe
    outfile = f"{cohort}_mergedMethylation_{datetime.now().strftime('%Y-%m-%d')}.gzip.tsv.gz"
    log_time(f'write out df to {outfile}')
    # Write output to  file
    combined_meth.to_csv(f'{output_dir}/{outfile}', sep="\t", index=False, header=True, chunksize=100000, compression="gzip")
    subprocess.run(f"mv {output_dir}/{outfile} {output_dir}/{outfile}.tmp.gz && bgzip -d {output_dir}/{outfile}.tmp.gz | bgzip > {output_dir}/{outfile}", shell=True, check=True)


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
        help="Path to the hap2 cohort combined input TSV file with header."
    )

    parser.add_argument(
        "-c","--input_ungrouped_cohort_bed",
        type=str,
        required=True,
        help="Path to the ungrouped cohort combined input TSV file with header."
    )

    parser.add_argument(
        "-o","--output_dir",
        type=str,
        required=True,
        help="output directory."
    )

    parser.add_argument(
        "-n","--cohort_name",
        type=str,
        required=True,
        help="output cohort name for naming output files."
    )

    # add optional argument for the chromosomes to run on
    parser.add_argument(
        "-l", '--chromosomes', 
        nargs='*', 
        default=None, 
        help="List of chromosomes to process (optional). If not provided, all chromosomes will be used.")

    # add argument for threads to use, number of chr's to run at once
    parser.add_argument(
        "-t","--threads",
        type=int,
        default=3,
        help="threads; default 3."
    )

    parser.add_argument(
        "-m","--min_cov",
        type=int,
        required=False,
        default=5,
        help="Minimum coverage; default=5x."
    )

    parser.add_argument(
        "-w","--write_croms",
        type=bool,
        required=False,
        default=False,
        help="Boolean to determine if filtered positions per chromosome are written out; default False"
    )


    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Parse arguments
    args = parser.parse_args()

    # Create output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    log_time_top(f"Output directory ensured at: {output_dir}")

    # if no chromosomes are provided, use all chromosomes from the BED
    if args.chromosomes is None:
        # get chromosome list from tabix 
        args.chromosomes = subprocess.check_output(["tabix", "-l", args.input_h1_cohort_bed], text=True).splitlines()

    log_time_top(f'running methylation merge with min cov: {args.min_cov}')

    process_chromosomes_in_parallel(output_dir, args.input_h1_cohort_bed, args.input_h2_cohort_bed, args.input_ungrouped_cohort_bed, args.cohort_name, int(args.min_cov), args.write_croms, args.chromosomes, args.threads)

    # combine_haplotypes(chrom, args.input_h1_cohort_bed, args.input_h2_cohort_bed, args.input_ungrouped_cohort_bed, args.output_file, int(args.min_cov) )





