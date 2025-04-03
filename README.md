# CARD_RegionalMethylation

Methods to aggregate methylation from ONT R9 or R10 sequencing data after merging coverage and average methylation data across a cohort. 
Merging is done as a part of the [cohort-level](https://github.com/nanoporegenomics/napu_wf/blob/R10_gt/wdl/cohort_wdl/scripts/merge_modkit_beds_allCpGs.py) work of the NAPU pipeline. 

## Python3 Dependancies
 ```argparse, numpy, pandas, datetime, matplotlib, seaborn, sklearn, statsmodels```  

## Methylation Aggregation
In development 

## Multivariate Regression with Age
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
### Test Regression 
```
python3.9 methylAgeRegresser.py -m data/sample_covs.10.csv -b data/sample_cgi.bed -o test_out -c test_CGIs --num_pcs 0
```
There are no significant hits in the provided test data. It should produce 2 BED files, 2 csv files, and 4 images.

