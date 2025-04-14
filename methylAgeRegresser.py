#!/usr/bin/env python3

import pandas as pd
import os
import sys
import gzip
import gcsfs
from datetime import datetime
import argparse

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import math
from numpy.polynomial import Polynomial
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns


"""
Performs a regression Meth ~ age + sex + PMI

inputs: (1) a metadat/covariate file with at least Age, Gender, and PMI 
            This method has been developed for the NABEC and HBCC cohorts
        (2) a aggregated methylBED file with the first 3 columns being chrom, start, end folowed by samples as columns
            Some aggregated methylBEDs have extra columns, UCSC CGI's for example 

example: python3 methylAgeRegressor.py -i NABEC_cohort_methyl_012025.tsv -o out_HarmPhase_methylationBeds

Outputs: 
        (1) Scatter plot of PC1 & PC2
        (2) Barchart of variance explalined by all PCs
        (3) Output of multivariable linear regression as a .csv
        (4) Significant hits ( FDR < alpha ) of multivariable linear regression as a .csv
        (5) Significant hits ( FDR < alpha ) of multivariable linear regression as a BED file with Age coeff as value
        (6) Significant hits ( FDR < alpha ) of multivariable linear regression as a BED file with Age coeff as value
            for regions with abs(AgeCoeff)>0.1 and r2>0.4
        (7) Volcano plot of all significant hits..

        
Author: Melissa Meredith
4/2025
"""

def log_time(message):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}")

def get_meta(metafile):
    """ load in metadata for either cohort, 
        meta columns are different for each 
        reformat ids to match bed
    """

    # ensure the delimiter matches
    delimiter=","
    if metafile[-4:]==".tsv":
        delimiter="\t"

    covs = pd.read_csv(metafile, sep=delimiter)

    sample_Id_Col = list(covs.columns)[0]
    print('sample_Id_Col', sample_Id_Col) 

    # if the sample ID's don't end with FTX add the ending to match the methylation sample ids
    if not covs[sample_Id_Col].iloc[0].endswith("_FTX"):
        covs[sample_Id_Col] = covs[sample_Id_Col]+"_FTX"
    covs.set_index(sample_Id_Col, inplace=True)

    # if the age and gender columns aren't named Age and Gender change them
    col_names = list(covs.columns)
    if ['Gender', 'Age'] not in col_names:

        # This is really just for the HBCC cohort naming scheme
        if "AgeDeath" in col_names:
            covs.rename(columns={"Sex at Birth": "Gender", "AgeDeath":"Age"}, inplace=True)
            # drop these extra columns as well for HBCC
            covs.drop(["Ethnicity", "Race"], axis=1, inplace=True)
            # one-hot encode gender HBCC is integers already
            covs['Gender'] = covs['Gender'].str.lower()
            # covs = pd.get_dummies(covs, columns=['Gender'], drop_first=True)

    # else:
    # one-hot encode gender for the NABEC cov data 
    covs = pd.get_dummies(covs, columns=['Gender'], drop_first=True)  
    covs['Gender_male'] = covs['Gender_male'].astype(int)

    return covs


def colinearity(df, cohort_region, output_dir):
    """ check correlation between factors"""
    log_time(f"Correlations")

    coor_matrix = df.iloc[:,:100].corr()

    plt.figure(figsize=(18, 16))  
    sns.heatmap(coor_matrix, annot=False, fmt=".3f", cmap='coolwarm', vmin=-1, vmax=1, square=True)
    plt.title('Correlation Heatmap of Variables')
    plt.savefig(f"{output_dir}/{cohort_region}_correlationHeatmap_{datetime.now().strftime('%Y-%m')}.png", dpi=300)

def pca(methyl_autosomes, covs, cohort_region):
    """ runs a principal component analysis of the input bed methyl data """
    
    log_time(f"Run PCA (20 components), plot variance explained, and PC1 & PC2")
    print(methyl_autosomes.T.head())

    # standardize the data; mean=0, variance=1
    scaler = StandardScaler()
    # Run pca on each sample to account for variable methylation across samples?
    # x_scaled = scaler.fit_transform(methyl_autosomes.dropna(axis=0))

    # Run pca on each region to capture some latent variable in the data
    x_scaled = scaler.fit_transform(methyl_autosomes.T.dropna(axis=1))

    # Innicialize pca of methylation data
    n_components=20
    if n_components>x_scaled.shape[0]:
        n_components=x_scaled.shape[0]-1
    pca = PCA(n_components=n_components)

    # fit the pca to the data's variation     
    pca.fit(x_scaled)
    pca_row = pca.transform(x_scaled)


    # plot variance explained 
    fig, axs = plt.subplots(1,2, figsize=(15,6))
    axs[0].bar(range(1, n_components+1), pca.explained_variance_ratio_ )
    axs[0].set_xticks(range(1, n_components+1))
    axs[0].set_title("PCA variance explained")
    axs[0].set_ylabel("Variance")
    axs[0].set_xlabel("PC")

    cumulative_variance = np.cumsum(pca.explained_variance_ratio_)
    axs[1].plot(range(1, n_components+1), cumulative_variance, marker='o', linestyle='--')
    axs[1].axhline(y=0.95, color='r', linestyle='--', label='95% Explained Variance')
    axs[1].legend(loc='upper right')
    axs[1].set_xticks(range(1, n_components+1))
    axs[1].set_title("PCA cumulative variance explained")
    axs[1].set_ylabel("Variance")
    axs[1].set_xlabel("PC")

    plt.savefig(f"{output_dir}/{cohort_region}_meth_pca_varianceExplained_20components_{datetime.now().strftime('%Y-%m')}.png", dpi=300)

    # plot first 2 pca's
    pca_df = pd.DataFrame( data=pca_row[:,:2], columns=['PC1', 'PC2'])
    pca_df = pd.concat([pca_df, covs.reset_index()], axis=1)

    plt.figure(figsize=(8,8))
    sns.scatterplot(x="PC1", y="PC2",data=pca_df, hue="Age", size="Gender_male")
    plt.axis('equal') 
    plt.title("WG Meth PCA")
    plt.xlabel(f"PC1 : {round(pca.explained_variance_ratio_[0],4)} variance explained", fontsize=14)
    plt.ylabel(f"PC2 : {round(pca.explained_variance_ratio_[1],4)} variance explained", fontsize=14)

    plt.savefig(f"{output_dir}/{cohort_region}_meth_pca_20components_{datetime.now().strftime('%Y-%m')}.png", dpi=300)

    # name the PCs in a dataframe
    pc_names = [f'PC{i+1}' for i in range(n_components)]
    pca_df = pd.DataFrame(pca_row, columns=pc_names)

    return pca_df


def regression_covs(data_df, num_pcs, cohort_region, date):
    """ regression function including PCs as covariates that capture variance across samples and use  """
    log_time(f"Linear Regression")

    if num_pcs > 0:
        included_pc_columns = [f'PC{i+1}' for i in range(num_pcs)]
    else: 
        included_pc_columns = []

    region_lreg = {}

    for column in data_df.columns:
        if column not in ['Age','Gender_male','GiB_of_uBAM','PMI','Region']+included_pc_columns:
            clean_data_df = data_df[['Age', 'PMI', 'Gender_male', column]+included_pc_columns].dropna()

            X = clean_data_df[['Age', 'PMI', 'Gender_male']+included_pc_columns]  # I variables
            y = clean_data_df[[column]]  # D variable must be 2D array  
            # if there are no rows, skip it
            if X.shape[0] == 0:
                continue
        
            # add a constant term (intercept) to the independent variable
            X = sm.add_constant(X)
    

            # fit the regression model
            model = sm.OLS(np.asarray(y), np.asarray(X)).fit()

            # print(model.summary() )
            
            # get p-values  
            region_lreg[column] = { 'r2': model.rsquared,
                                      'intercept':model.params[0],
                                      'Age_Coeff': model.params[1], 
                                      'PMI_Coeff': model.params[2],
                                      'Gender_coeff':model.params[3],
                                      'P_Value_Age': model.pvalues[1],
                                      'P_Value_PMI': model.pvalues[2],
                                      'P_Value_Gender': model.pvalues[3],
                                      
                                     }

    # store the LR output as dataframe, separate parameters, and write to csv/bed
    data_cov_lr_df = pd.DataFrame(region_lreg).T
    
    data_cov_lr_df.to_csv(f'{output_dir}/{cohort_region}_olsRcov_Age_PMI_Gender_{date}.csv', sep='\t', header=True,index=True)

    return data_cov_lr_df

def plot_volcano(df, cohort_region, alpha):
    """ Linear Regression Volcano Plot""" 

    log_time(f"Plotting volcano")

    fig, axs = plt.subplots(1, 1, figsize=(8, 10)) 

    sns.scatterplot(x='Age_Coeff', y='neg_log10_bh_corrected_P_Value_Age', data=df, hue='neg_log10_bh_corrected_P_Value_Age', edgecolor='k', s=28, ax=axs, palette='coolwarm')
    bh_num_sig = df.loc[(df['bh_corrected_P_Value_Age']<0.05)].shape[0]
    axs.axhline(y=-np.log10(alpha), color='red', label=f'FDR=0.05: {bh_num_sig} sig hits')

    axs.legend()
    axs.set_title(f'{cohort_region} Avg Methylation Linear Regression x Age')
    axs.set_ylabel("B-H Corrected - log10(p)")
    axs.set_xlabel("Slope of Avg Methylation x Age Linear Regression")
    axs.grid(True)

    # t = plt.xticks(rotation=90)

    plt.savefig(f"{output_dir}/{cohort_region}_Age_r2_x_LinearRegressionSlopeVolcano_Benjamini-Hochberg_{datetime.now().strftime('%Y-%m')}.png",dpi=300, facecolor='white', transparent=False)



def format_run_regression(methylBed, covs, haplotype, alpha, num_pcs, cohort_region, output_dir):
    """ set up data for regression and run it """

    log_time(f"Read in methylbed")
    # load in the methyl bed and replace . with na values 
    meth_bed = pd.read_csv(methylBed,  delimiter='\t', na_values=['.'])

    # make a position column out of chrom, start, and end
    # check for chrom syntax
    chrom_name = "chrom"
    if "#chrom" in meth_bed.columns:
        chrom_name = "#chrom"

    start_col = "start"
    if start_col not in meth_bed.columns:
        start_col = meth_bed.columns[1]

    end_col = "end"
    if end_col not in meth_bed.columns:
        end_col = meth_bed.columns[2]

    # make a position column and set it as the index
    meth_bed['position'] = meth_bed[chrom_name] + '_' + meth_bed[start_col].astype(str) + '_' + meth_bed[end_col].astype(str)
    meth_bed.drop([chrom_name, start_col,end_col], axis=1, inplace=True)
    meth_bed.set_index(['position'], inplace=True)
    # remove any extra columns
    # meth_bed = meth_bed.loc[:, meth_bed.columns.str.startswith('avgMod')]
    meth_bed = meth_bed.loc[:, meth_bed.columns.str.endswith('_modFraction')]

    # replace avgMod_sampleID with just sampleIDs to match covariates 
    # meth_bed.rename(columns=lambda x: x.replace('avgMod_', '') if isinstance(x, str) and x.startswith('avgMod_') else x, inplace=True)
    # check which haplotype and replace column names to match covariates
    str_to_replace = ""
    if haplotype == 1:
        str_to_replace = '_GRCh38_1_modFraction'
    elif haplotype == 2:
        str_to_replace = '_GRCh38_2_modFraction'
    else:
        log_time(f"Haplotype {haplotype} not 1 or 2?; Exiting")
        return -1

    meth_bed.rename(columns=lambda x: x.replace(str_to_replace, '') if isinstance(x, str) and x.endswith('_modFraction') else x, inplace=True)
    print( meth_bed.head())

    log_time(f"Drop non autosomes")
    # drop sex chromosomes
    meth_bed_autosomes = meth_bed[~meth_bed.index.str.contains('X|Y')]
    # ensure all values are floats and drop positions with a NA for any sample
    meth_bed_autosomes = meth_bed_autosomes.astype(float)
    meth_bed_autosomes = meth_bed_autosomes.dropna(axis=0)
    log_time(f"Total Regions:{meth_bed.shape[0]}. {meth_bed_autosomes.shape[0]} regions after dropping regions with NaN values")

    # Run pca
    pca_df = pca(meth_bed_autosomes, covs, cohort_region)

    log_time(f"Join methyl with covarites and {num_pcs} PCs")
    # join methylation and covariates
    meth_bed_autosomes_age = covs.join(meth_bed_autosomes.T, how = 'inner')
    meth_bed_autosomes_age_pcs = pd.concat([pca_df.iloc[:,:num_pcs], meth_bed_autosomes_age.reset_index()], axis=1)
    meth_bed_autosomes_age_pcs.set_index('index', inplace=True)
    print(meth_bed_autosomes_age_pcs.head())
    colinearity(meth_bed_autosomes_age_pcs.drop(['Region'], axis=1), cohort_region, output_dir)



    # Perform linear regressions
    meth_bed_autosomes_age_lr_df = regression_covs(meth_bed_autosomes_age_pcs, num_pcs, cohort_region, datetime.now().strftime('%Y-%m') )
    print(meth_bed_autosomes_age_lr_df.head())
    
    log_time(f"B-H Correction")
    # Apply Benjamini-Hochberg correction
    reject, p_values_corrected, _, fdr = multipletests(meth_bed_autosomes_age_lr_df['P_Value_Age'], alpha=alpha, method='fdr_bh')


    meth_bed_autosomes_age_lr_df['bh_corrected_P_Value_Age'] = p_values_corrected
    meth_bed_autosomes_age_lr_df['neg_log10_bh_corrected_P_Value_Age'] = -np.log10(p_values_corrected)

    log_time(f"Write out significant hits") 
    # Reformat the chr position to bed format 
    meth_bed_autosomes_age_lr_df_noi = meth_bed_autosomes_age_lr_df.reset_index()
    meth_bed_autosomes_age_lr_df_noi.rename(columns={"index": "Position"}, inplace=True)
    print(meth_bed_autosomes_age_lr_df_noi.head())

    meth_bed_autosomes_age_lr_df_noi.loc[meth_bed_autosomes_age_lr_df_noi['bh_corrected_P_Value_Age']<alpha].to_csv(f"{output_dir}/{cohort_region}_olsRcov_Age_PMI_Gender_BHsig_{datetime.now().strftime('%Y-%m')}.csv",index=True,header=True,sep=",")

    df_split = meth_bed_autosomes_age_lr_df_noi['Position'].str.split('_', expand=True)
    df_split.columns = ['#chrom', 'start','end']

    meth_bed_autosomes_age_lr_df_noi['#chrom'] = df_split['#chrom']
    meth_bed_autosomes_age_lr_df_noi['start'] = df_split['start']
    meth_bed_autosomes_age_lr_df_noi['end'] = df_split['end']


    # write out all the B-H Significant hits; and then those with slope (-.05 > S > 0.05)
    meth_bed_autosomes_age_lr_df_noi_bhSig = meth_bed_autosomes_age_lr_df_noi.loc[meth_bed_autosomes_age_lr_df_noi['bh_corrected_P_Value_Age'] < 0.05]

    log_time(f"Write out all hits with bh- corrected p value as bed")
    meth_bed_autosomes_age_lr_df_noi_bhSig[['#chrom', 'start','end','bh_corrected_P_Value_Age','Age_Coeff']].to_csv(f'{output_dir}/{cohort_region}_lr_Age_PMI_Gender.bed',index=False,header=False,sep="\t")
    log_time(f"Write out all hits with bh- corrected slope as bed")
    meth_bed_autosomes_age_lr_df_noi_bhSig.loc[(meth_bed_autosomes_age_lr_df_noi_bhSig['Age_Coeff'] > 0.1) & (meth_bed_autosomes_age_lr_df_noi_bhSig['r2'] > 0.4)][['#chrom', 'start','end','Age_Coeff']].to_csv(f'{output_dir}/{cohort_region}_lr_Age_PMI_Gender_bh_sig_AgeSlope1_r2.4.bed',index=False,header=False,sep="\t")

    plot_volcano(meth_bed_autosomes_age_lr_df_noi_bhSig, cohort_region, alpha)






if __name__ == "__main__":

    # Make an argument parser
    parser = argparse.ArgumentParser(description="Merge methylBeds from Modkit across cohorts.")
    parser.add_argument(
        "-m","--in_meta_file",
        type=str,
        required=True,
        help="Path to the input metadata/covariate file."
    )

    parser.add_argument(
        "-o","--output_directory",
        type=str,
        required=True,
        help="name of output directory."
    )

    parser.add_argument(
        "-b","--methylbed",
        type=str,
        required=True,
        help="methylbed file with samples as columns"
    )

    parser.add_argument(
        "-c","--cohort_region",
        type=str,
        required=True,
        help="cohort name and region name for naming output files"
    )

    parser.add_argument(
        "-d","--haplotype",
        type=int,
        required=True,
        help="haplotype; integer 1 or 2"
    )

    parser.add_argument(
        "-a","--alpha",
        type=float,
        default=0.05,
        required=False,
        help="alpha significance value for FDR, default:0.05"
    )

    parser.add_argument(
        "-p","--num_pcs",
        type=int,
        default=10,
        required=False,
        help="Number of principal components to include in regression, default:10"
    )


    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Parse arguments
    args = parser.parse_args()

    
    # Create output directory
    output_dir = args.output_directory
    os.makedirs(output_dir, exist_ok=True)
    log_time(f"Output directory ensured at: {output_dir}")

    # Read in covariates 
    cov_data = get_meta(args.in_meta_file)

    # Set up and run regression 
    format_run_regression(args.methylbed, cov_data, args.haplotype, args.alpha, args.num_pcs, args.cohort_region, output_dir)





