# CARD_RegionalMethylation

Methods to aggregate methylation from ONT R9 or R10 sequencing data after merging coverage and average methylation data across a cohort. <br>
Merging is done as a part of the [cohort-level](https://github.com/nanoporegenomics/napu_wf/blob/R10_gt/wdl/cohort_wdl/scripts/merge_modkit_beds_allCpGs.py) work of the NAPU pipeline. 

## Dependancies
``` bedtools, bgzip, tabix, python3```
## Python3.9 Dependancies
 ```argparse, numpy, pandas, datetime, matplotlib, seaborn, sklearn, statsmodels```  
 

## Methylation Aggregation
Methylation is averaged across a given set of BED regions of interest. <br>
The process is automated in the ``` organizeSortedMethylBed.sh ``` , but the individual steps are also described below. 

### Filtering and Aggregation Bash Script
The input cohort BED/tsv and input regional BED need to be sorted, bgzip-ed, and indexed for parallel analysis. <br>
The ``` organizeSortedMethylBed.sh ``` bash script automates these steps of sorting, bgzipping, indexing. It also runs the filtering script and the bedtools map aggregation step. 

Change permissions to be executable: <br>
```chmod u+x organizeSortedMethylBed.sh ``` <br>
Test executable: <br>
```
./organizeSortedMethylBed.sh
```
Expected output:
```

 Usage: ./organizeSortedMethylBed.sh <region.bed> <cohort.tsv> <filter.py> <cohortId> 
```
filter.py is the location of the ```filter_cohort_methylation_parallel.py``` script. <br>
Each of the commands in the bash script could be run individually if desired. 

Example Usage with CpG Islands BED: <br>
``` 
$ ./organizeMethylBed.sh small.cpg.bed cohort.small.tsv.gz filter_cohort_methylation_parallel.py <COHORT_NAME>

```
<br>

The whole genome CpGs are intersected with the bed regions of interst in order to reduce the number of positions analyzed downstream.  <br>

### Filtering and Aggregation Bash Script Output

```
(1) input group BED intersected with bed region: cohort.small.tsv.gz.small.cpg.bed.tsv.gz  
(2) intersected index:                           cohort.small.tsv.gz.small.cpg.bed.tsv.gz.tbi  
(3) bgzip-ped and indexed bed regions:           small.cpg.bed.sorted.tsv.gz  small.cpg.bed.sorted.tsv.gz.tbi  
(4) directory for coverage/cpg filtered output:  filter_small.cpg.bed_region   
(5) File containing the header for mapping:      regional_bedtoolsMapMean_header.tsv
(6) File with mappings for each region:          cohort.small.tsv.gz.small.cpg.bed.regional_bedtoolsMapMean.cohort.tsv
```
The file with mappings for each region ```cohort.small.tsv.gz.small.cpg.bed.regional_bedtoolsMapMean.cohort.tsv```is what can be used for regression. <br>
<br>
The following sections describe the filtering that is performed inside of the bash script:

#### Filtering and Aggregation Inputs 

CpG positions that intersect with the input regional BED are considered for minimum coverage and minimum number of CpGs within that region. 
```filter_cohort_methylation_parallel.py``` removes individual CpGs for entire regions where these minimums are not covered. <br> <br>
Filtering individual CpG sites that dropout is under development, currently an average coverage is used to evalulate passing coverage minimum. 

```
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
```
There are 2 required inputs: A cohort BED/tsv and a regional BED file. <br><br>
The cohort BED/tsv input is a phased BED/tsv file for each CpG in the reference genome (GRCh38)<br>
including the usual BED positions ```#chrom start end``` <br>
followed by 3 columns for each sample in a cohort: <br>
```SAMPLE_GRCh38_2_validCov        SAMPLE_GRCh38_2_modFraction     SAMPLE_GRCh38_2_modReads  ```

The second input required is a BED file of regions of interest. <br>
For example CpG Islands in GRCh38: 
```
data/small.cpg.bed 
#chrom  chromStart  chromEnd  name  length  cpgNum  gcNum perCpg  perGc obsExp
chr1  28735 29737 CpG:_111  1002  111 731 22.2  73  0.85
chr1  135124  135563  CpG:_30 439 30  295 13.7  67.2  0.64
chr1  199251  200121  CpG:_104  870 104 643 23.9  73.9  0.89
```

#### Filtering outputs:
A tsv of average coverage per region, positions that failed filters, and a file of positions that passed filters. <br> 

```
cohort.small_hap1_tsv_cpg_hg38_bed_tsv_average_regional_coverage_2025-04-12.tsv  
cohort.small_hap1_tsv_cpg_hg38_bed_tsv_failedCovFilter_5cov_5cpgs_2025-04-12.tsv  
cohort.small_hap1_tsv_cpg_hg38_bed_tsv_filtered_5cov_5cpgs_2025-04-12.tsv
```
<br> 
The file of positions that passed filters is then passed to bedtools map to average methylation over the input bed regions. 

```
bedtools map -a small.cpg.bed.sorted.tsv.gz -b filter_small.cpg.bed_region/cohort_small_tsv_gz_small_cpg_bed_tsv_filtered_5cov_5cpgs_2025-04-15.tsv -c  5,8,11 -o mean 
```

The output of the bedtools map command can then be used for age regressions. <br>


## Multivariate Regression with Age
After filtering for coverage and cpg count the aggregated methylation file can be used for regression  ( the "regional_bedtoolsMapMean" file ). 
``` methylAgeRegresser.py``` Runs a multivariate regression on aggregated methylation data on all samples in a cohort. 
The regression equation is: 
``` Methylation ~ Age + Gender + PostMortemInterval + PCs (optional) ``` 

### Regression Inputs
```
Required:
  -m IN_META_FILE, --in_meta_file IN_META_FILE
                        Path to the input metadata/covariate file.
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        name of output directory.
  -b METHYLBED, --methylbed METHYLBED
                        methylbed file with samples as columns
  -c COHORT_REGION, --cohort_region COHORT_REGION
                        cohort name and region name for naming output files
  -d HAPLOTYPE, --haplotype HAPLOTYPE
                        haplotype; integer 1 or 2; or 0 for unphased
Options:
  -a ALPHA, --alpha ALPHA
                        alpha significance value for FDR, default:0.05
  -p NUM_PCS, --num_pcs NUM_PCS
                        Number of principal components to include in regression; default 10
```

### Regression Outputs
```
(1) Scatter plot of PC1 & PC2
(2) Barchart of variance explalined by all PCs
(3) Heatmap of correlations for covariates and first ~50 features
(4) Output of multivariable linear regression as a .csv
(5) Significant hits ( FDR < alpha ) of multivariable linear regression as a .csv
(6) Significant hits ( FDR < alpha ) of multivariable linear regression as a BED file with Age coeff as value
(7) Significant hits ( FDR < alpha ) of multivariable linear regression as a BED file with Age coeff as value
  for regions with abs(AgeCoeff)>0.1 and r2>0.4
(8) Volcano plot of all significant hits.
```
### Test Phased Regression 
```
python methylAgeRegresser.py -m data/sample_covs.10.csv -b data/small.phased.testData.tsv -o test_phased -d 2 -c test_data --num_pcs 0
```
There are no significant hits in the provided test data. It should produce 2 BED files, 2 csv files, and 4 images.
Example output is in ```data/test_phased/```

### Test Unphased Regression ( depreciated )
```
python methylAgeRegresser.py -m data/sample_covs.10.csv -b data/sample_cgi.bed -o test_out -d 0 -c test_CGIs --num_pcs 0
```
There are no significant hits in the provided test data. It should produce 2 BED files, 2 csv files, and 4 images.
Example output is in ```data/test_out/```

