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
sort
(head -n 1 $cohorttsv && tail -n +2 $cohorttsv | sort -k 1,1 -k2,2n ) > $clean_intsv_name.sorted.tsv



    Filter a cohort modkit tsv by valid coverage, for aggregating over bed regions. Input bed and cohort tsv should
    be bgzip and tabix'ed, the index is used for quick lookup when running chromsomes in parallel. 
    Also it is faster to intersect the cohort tsv/bed with the region bed for less regions to look at: 
    bedtools intersect -header -wa -a combined_methylation_hap1.tsv.sorted.tsv -b ../small.cpg.bed > indiv.cpgs.inregion.tsv


    usage: filter_cohort_methylation_parallel.py [-h] -i INPUT_COHORT_TSV -b IN_BED_FILE -o OUTPUT_DIR

    Filter positions from the group BED that don't pass coverage or cpg minimum filters.

    Required arguments:
      -h, --help            show this help message and exit
      -i INPUT_COHORT_TSV, --input_cohort_tsv INPUT_COHORT_TSV
                            Path to the cohort combined input TSV file with header ( must be bgzip and tabix ).
      -b IN_BED_FILE, --in_bed_file IN_BED_FILE
                            Path to the bed file region of interest ( must be bgzip and tabix ).
      -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                            output directory.

    Optional arguments:
      -l [CHROMOSOMES ...], --chromosomes [CHROMOSOMES ...]
                            List of chromosomes to process (optional). If not provided, all chromosomes will be used.
      -m MIN_COV, --min_cov MIN_COV
                            Minimum coverage; default=5x.
      -g MIN_CPGS, --min_cpgs MIN_CPGS
                            Minimum number of CpGs within a bed region; default=5x.
      -c AVERAGE_COVERAGE, --average_coverage AVERAGE_COVERAGE
                            Boolean to determine if every cpg site needs to be > min cov or average across region is > mincov; default True
      -f WRITE_FAILS, --write_fails WRITE_FAILS
                            Boolean to determine if the positions that don't pass coverage or min cpg filter are written out; default True
      -t THREADS, --threads THREADS
                            threads; default 3.

"""


def log_time(message):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}")

def check_min_cov(val_list, min_cov):
    """ check if all coverage means are at >= the min_cov cutoff, 
        if any don't return the indices of those values.
    """
    # log_time(f'checking for coverage')

    low_cov_indices = [i for i,val in enumerate(val_list) if val < min_cov]

    # if values didn't pass cov min return false and list of indices, otherwise return True
    if low_cov_indices: 
        return False, low_cov_indices
    else:
        return True, []

def mask_mod_fraction(row, validCov_cols, modFraction_cols, min_cov):
    for cov_col, mod_col in zip(validCov_cols, modFraction_cols):
        if row[cov_col] < min_cov:
            row[mod_col] = np.nan
    return row

def filter_regions_by_coverage(chrom, in_tsv, in_bed, average_coverage, min_cov, min_cpgs):
    """ average coverage for each input bed region and store pass and failed regions
        
        Consider a method to remove samples with poor coverage..

        subset region by chromosome being run here
    """

    log_time(f'process {chrom}')
    # query the chromosome from the bgzipped BED file
    cmd = ["tabix", in_tsv, chrom]
    try:
        # log_time(f'command {cmd}')
        in_tsv_chrom = subprocess.check_output(cmd, text=True)
        # log_time(f'zcat {in_tsv}')
        header = subprocess.check_output(f"zcat {in_tsv} | head -1", shell=True, text=True).strip().split("\t")     ## fix before leaving local 
        # header = subprocess.check_output(f"gzip -dc {in_tsv} | head -1", shell=True, text=True).strip().split("\t")
    except subprocess.CalledProcessError as e:
        print(f"Failed tsv on {chrom}: {e}")
        return chrom, 0

    try:
        log_time(f'load in {chrom} subset')
        rdata = pd.read_csv(StringIO(in_tsv_chrom), sep='\t')
        rdata.columns = header
        # if in_tsv[-3:] == ".gz":
        #     with gzip.open(in_tsv,'rb') as file1:
        #         rdata = pd.read_csv(file1, sep='\t')
        # else:
        #     rdata = pd.read_csv(in_tsv, sep='\t')
    except Exception as e:
        log_time(f"Error reading {chrom} files {e}" )
        return [], [], []

    # store all passsing regions in a list
    pass_regions = [] 
    fail_regions = []
    region_avg_coverages = []

    # store passing region averages
    # pass_region_averages = []




    log_time(f'filter for coverages')

    # Get column names that end with '_validCov'
    validCov_cols = [col for col in rdata.columns if col.endswith("_validCov")] 
    log_time(f'{len(validCov_cols)} samples')
    # get mod fraction columns 
    modFraction_cols = [col for col in rdata.columns if col.endswith("_modFraction")] 

    log_time(f'get {chrom} from {in_bed}')
    # query the chromosome from the bgzipped BED file
    cmd = ["tabix", in_bed, chrom]
    try:
        in_bed_chrom = subprocess.check_output(cmd, text=True)
        bedhead = subprocess.check_output(f"zcat {in_bed} | head -1", shell=True, text=True).strip().split("\t")       ## fix before leaving local 
        # bedhead = subprocess.check_output(f"gzip -dc {in_bed} | head -1", shell=True, text=True).strip().split("\t")
    except subprocess.CalledProcessError as e:
        print(f"Failed bed on {chrom}: {e}")
        return chrom, 0

    # write out the header for bedtools map:
    if (chrom=="chr1"):
        modFraction_cols = [col for col in rdata.columns if col.endswith("_modFraction")] 
        average_header = ("\t").join(bedhead) +"\t"+ ("\t").join(modFraction_cols) + "\n"
        with open("regional_bedtoolsMapMean_header.tsv", 'w') as outf:
                    outf.write(average_header)

    # with open(in_bed) as b:
    for line in in_bed_chrom.strip().split("\n"):
        if line.strip() == "" or line.startswith("#"):
            continue

        # get bed regions
        tokens = line.strip().split()[:3]
        chrom = tokens[0]
        start = int(tokens[1])
        end = int(tokens[2]) 

        # log_time(f'filtering {chrom}: {start} {end}')

        # get just the cpgs within the bed region 
        # TODO make sure this is #chrom?
        region_rdata = rdata.loc[ (rdata['#chrom']==chrom) &
                                  (rdata['start'] >= start) &
                                    (rdata['end'] <= end) ]

        
        # log_time(f'subset region. df shape: {region_rdata.shape}')
        # check Number of cpgs passes filter
        if region_rdata.shape[0] > min_cpgs:
            # print('region_rdata',start, end, region_rdata.iloc[:, :6])

            # if average is selected find avg coverage, default is false
            if average_coverage:

                log_time(f'find avg coverages')
                # make a list to store all the coverage averages per sample 
                region_coverage_averages = []
                for sample_cov_col in validCov_cols:
                    # get average coverage for the region per sample, round to 2 decimal points
                    # cov_mean = region_rdata[sample_cov_col].mean()
                    region_coverage_averages.append(np.round(region_rdata[sample_cov_col].mean(),2))

                # check if coverage mean passes fiilter
                pass_bool, failing_indices = check_min_cov(region_coverage_averages, min_cov)

                region_avg_coverages.append( pd.DataFrame([region_coverage_averages], columns=validCov_cols, index=[f'{chrom}_{start}_{end}']) )

                # if all pass store region in pass list
                if pass_bool:
                    # store all pass reigions 
                    pass_regions.append(region_rdata)
                else:
                    fail_regions.append(region_rdata)

                

            else:

                # remove rows that don't have valid cov for all samples
                # region_rdata_pass = region_rdata[
                    #     region_rdata.apply(lambda row: all(row[col] >= min_cov for col in validCov_cols), axis=1)
                    # ]
                # replace methylation values with np.nan if cov col is below threshold
                # region_rdata[validCov_cols] = region_rdata[validCov_cols].applymap(
                #                                 lambda x: x if x >= min_cov else np.nan
                #                                 )

                # region_rdata_pass = region_rdata.apply(
                #     lambda row: mask_mod_fraction(row, validCov_cols, modFraction_cols, min_cov),
                #     axis=1)

                # mask individual samples that don't have sufficient coverage 
                for cov_col, mod_col in zip(validCov_cols, modFraction_cols):
                    mask = region_rdata[cov_col] < min_cov
                    region_rdata.loc[mask, mod_col] = "." #np.nan
                region_rdata_pass = region_rdata
                
                # store the cpg sites that 
                region_rdata_fail = region_rdata[
                        region_rdata.apply(lambda row: all(row[col] < min_cov for col in validCov_cols), axis=1)
                    ]

                # might not be necessary 
                region_coverage_averages = []
                for sample_cov_col in validCov_cols:
                    # get average coverage for the region per sample
                    # cov_mean = region_rdata[sample_cov_col].mean()
                    region_coverage_averages.append(np.round(region_rdata[sample_cov_col].mean(),2))
                region_avg_coverages.append( pd.DataFrame([region_coverage_averages], columns=validCov_cols, index=[f'{chrom}_{start}_{end}']) )
                
                # sore passing rows
                # pass_regions = pd.concat([pass_regions, region_rdata], ignore_index=True)
                if not region_rdata_pass.empty:
                    pass_regions.append(region_rdata_pass)

                if not region_rdata_fail.empty:
                    fail_regions.append(region_rdata_fail)


        else:
            fail_regions.append(region_rdata)

    # combine pass regions into df
    if len(pass_regions)>0:
        combined_pass_regions = pd.concat(pass_regions, ignore_index=True)
    else: 
        combined_pass_regions = pd.DataFrame()
    # combine fail regions 
    if len(fail_regions)>0:
        # log_time(f'combine failed filter df')
        combined_fail_regions = pd.concat(fail_regions, ignore_index=True)
    else: 
        combined_fail_regions = pd.DataFrame()

    # combine avg coverages
    region_avg_cov_df = pd.concat(region_avg_coverages, ignore_index=True)
    region_avg_cov_df.index.name = "#RegionID"


    return combined_pass_regions, combined_fail_regions, region_avg_cov_df



# process each chromosome in parallel
def process_chromosomes_in_parallel(output_dir, in_tsv, in_bed, average_coverage, min_cov, min_cpgs, chromosomes, nprocesses, write_fails):

    # create a pool of worker processes
    with Pool(processes=nprocesses) as pool:
        results = pool.starmap(filter_regions_by_coverage, [(chrom, in_tsv, in_bed, average_coverage, min_cov, min_cpgs) for chrom in chromosomes])

    # combine results
    log_time(f'process parallel results')

    # for i,result in enumerate(results):
    #     print(i,'pass\n', result[0])
    #     print(i,'fail\n', result[1])
    #     print(i,'avg\n', result[2])


    # process parallel results 
    # pass regions are the first df in all lists
    pass_regions_combined = pd.concat([result[0] for result in results], ignore_index=True)
    # fail regions are the second df in all lists
    fail_regions_combined = pd.concat([result[1] for result in results], ignore_index=True)
    # cov averages are the first df in all lists
    region_avg_cov_combined = pd.concat([result[2] for result in results], ignore_index=True)

    log_time(f'write out filtered df')
    # write out the filtered dataframe
    suffix = in_tsv.split("/")[-1]
    fileprefix = ("_").join(suffix.split(".")[:-1])

    if not pass_regions_combined.empty:
        # combined_pass_regions = pd.concat(pass_regions, ignore_index=True)
        outfile = f"{output_dir}/{fileprefix}_filtered_{min_cov}cov_{min_cpgs}cpgs_{datetime.now().strftime('%Y-%m-%d')}.tsv"
        pass_regions_combined.to_csv(outfile, sep="\t", index=False, header=True)
        print(f"Filtered regions written to {outfile} with {pass_regions_combined.shape[0]} rows.")
    else:
        print("No regions passed the filter.")

    if not fail_regions_combined.empty and write_fails:
        log_time(f'write out failed filter df')
        # comibined_fail_regions = pd.concat(fail_regions, ignore_index=True)
        outfile = f"{output_dir}/{fileprefix}_failedCovFilter_{min_cov}cov_{min_cpgs}cpgs_{datetime.now().strftime('%Y-%m-%d')}.tsv"
        fail_regions_combined.to_csv(outfile, sep="\t", index=False, header=True)

    log_time(f'write out average coverages per samples')
    # region_avg_cov_df = pd.concat(region_avg_coverages)
    region_avg_cov_combined.index.name = "#RegionID"
    region_avg_cov_combined.to_csv(f"{output_dir}/{fileprefix}_average_regional_coverage_{datetime.now().strftime('%Y-%m-%d')}.tsv", sep="\t", index=True, header=True)


if __name__ == "__main__":

    # Make an argument parser
    parser = argparse.ArgumentParser(description="Filter positions from the group BED that don't pass coverage or cpg minimum filters.")
    parser.add_argument(
        "-i","--input_cohort_tsv",
        type=str,
        required=True,
        help="Path to the cohort combined input TSV file with header ( must be bgzip and tabix )."
    )

    parser.add_argument(
        "-b","--in_bed_file",
        type=str,
        required=True,
        help="Path to the bed file region of interest ( must be bgzip and tabix )."
    )

    parser.add_argument(
        "-o","--output_dir",
        type=str,
        required=True,
        help="output directory."
    )

    # add optional argument for the chromosomes to run on
    parser.add_argument(
        "-l", '--chromosomes', 
        nargs='*', 
        default=None, 
        help="List of chromosomes to process (optional). If not provided, all chromosomes will be used.")

    parser.add_argument(
        "-m","--min_cov",
        type=int,
        required=False,
        default=3,
        help="Minimum coverage; default=3x."
    )

    parser.add_argument(
        "-g","--min_cpgs",
        type=int,
        required=False,
        default=5,
        help="Minimum number of CpGs within a bed region; default=3x."
    )

    parser.add_argument(
        "-c","--average_coverage",
        type=bool,
        required=False,
        default=False,
        help="Boolean to determine if every cpg site needs to be > min cov or average across region is > mincov; default False"
    )

    parser.add_argument(
        "-f","--write_fails",
        type=bool,
        required=False,
        default=True,
        help="Boolean to determine if the positions that don't pass coverage or min cpg filter are written out; default True"
    )

    # add argument for threads to use, number of chr's to run at once
    parser.add_argument(
        "-t","--threads",
        type=int,
        default=3,
        help="threads; default 3."
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Parse arguments
    args = parser.parse_args()

    # Create output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    log_time(f"Output directory ensured at: {output_dir}")

    # if no chromosomes are provided, use all chromosomes from the BED
    if args.chromosomes is None:
        # get chromosome list from tabix 
        args.chromosomes = subprocess.check_output(["tabix", "-l", args.in_bed_file], text=True).splitlines()

    log_time(f'running coverage filtering with: min cpgs: {args.min_cpgs} and min cov: {args.min_cov}')
    # filter_regions_by_coverage('chr1', args.input_cohort_tsv, args.in_bed_file, args.average_coverage, args.min_cov, args.min_cpgs)

    process_chromosomes_in_parallel(output_dir, args.input_cohort_tsv, args.in_bed_file, args.average_coverage, args.min_cov, args.min_cpgs, args.chromosomes, args.threads, args.write_fails)








