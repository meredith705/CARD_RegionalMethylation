#!/bin/bash

if [ "$#" -ne 4 ]; then
	echo "Usage: $0 <region.bed> <cohort.tsv> <filter.py> <cohortId>"
	exit 1
fi

bedregions=$1
cohorttsv=$2
filterScript=$3
cohort=$4

clean_intsv_name=$(basename "$cohorttsv")
clean_bed_name=$(basename "$bedregions")

# if [ ! -f "$clean_intsv_name.sorted.tsv.gz" ]; then
# 	# sort the input file
# 	echo "[$(date +"%Y-%m-%d %H:%M:%S")] Starting sort"
# 	(head -n 1 "$cohorttsv" && tail -n +2 "$cohorttsv" | sort -k 1,1 -k2,2n ) | bgzip -@16 > $clean_intsv_name.sorted.tsv.gz
# 	echo "[$(date +"%Y-%m-%d %H:%M:%S")] sorted and bgzip-ed"
# else
# 	echo "$clean_intsv_name.sorted.tsv.gz already exists."
# fi

if [ ! -f "$clean_intsv_name.$clean_bed_name.tsv.gz" ]; then
	echo "[$(date +"%Y-%m-%d %H:%M:%S")] intersect $cohorttsv with $bedregions "
	bedtools intersect -header -wa -a $cohorttsv -b $bedregions | bgzip -@16 > $clean_intsv_name.$clean_bed_name.tsv.gz
else
	echo "$clean_intsv_name.$clean_bed_name.tsv.gz already intersected "
fi

if [ ! -f "$clean_intsv_name.$clean_bed_name.tsv.gz.tbi" ]; then
	echo "[$(date +"%Y-%m-%d %H:%M:%S")] Starting index"
	tabix -p bed $clean_intsv_name.$clean_bed_name.tsv.gz
else
	echo "$clean_intsv_name.$clean_bed_name.tsv.gz.tbi already indexed "
fi

if [ ! -f "$clean_bed_name.sorted.tsv.gz.tbi" ]; then
	echo "[$(date +"%Y-%m-%d %H:%M:%S")] Sort, bgZip and index bed regions"
	(head -n 1 $bedregions && tail -n +2 $bedregions | sort -k 1,1 -k2,2n ) | bgzip -@8 > $clean_bed_name.sorted.tsv.gz
	tabix -p bed $clean_bed_name.sorted.tsv.gz
else
	echo "$clean_bed_name.sorted.tsv.gz.tbi already sorted and indexed "
fi

# Filter regions with low coverage in any sample
echo "[$(date +"%Y-%m-%d %H:%M:%S")] Starting coverage filtering"
outdirname="filter_${clean_bed_name}_region"
python $filterScript -o $outdirname -i ${clean_intsv_name}.${clean_bed_name}.tsv.gz -b ${clean_bed_name}.sorted.tsv.gz -t 12
cp regional_bedtoolsMapMean_header.tsv "$clean_intsv_name.$clean_bed_name.regional_bedtoolsMapMean.$cohort.tsv"

# look for output files
if [ -d "$outdirname" ]; then
	filtfile=$(ls "$outdirname"/*_filtered_* 2>/dev/null | head -n 1)

	if [ -z "$filtfile" ]; then
		echo "No file matching *_filtered_* found in $outdirname."
		exit 1
	else
		echo "Using file: $filtfile"
		bedtools sort -header -i ${filtfile} > ${filtfile}.sorted.tsv  #TEST THIS
		# check which cohort to run mapping 
		if [ "$cohort" == "NABEC" ]; then
			echo "[$(date +"%Y-%m-%d %H:%M:%S")] bedtools map NABEC"

			bedtools map -a ${clean_bed_name}.sorted.tsv.gz -b ${filtfile}.sorted.tsv -c 5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95,98,101,104,107,110,113,116,119,122,125,128,131,134,137,140,143,146,149,152,155,158,161,164,167,170,173,176,179,182,185,188,191,194,197,200,203,206,209,212,215,218,221,224,227,230,233,236,239,242,245,248,251,254,257,260,263,266,269,272,275,278,281,284,287,290,293,296,299,302,305,308,311,314,317,320,323,326,329,332,335,338,341,344,347,350,353,356,359,362,365,368,371,374,377,380,383,386,389,392,395,398,401,404,407,410,413,416,419,422,425,428,431,434,437,440,443,446,449,452,455,458,461,464,467,470,473,476,479,482,485,488,491,494,497,500,503,506,509,512,515,518,521,524,527,530,533,536,539,542,545,548,551,554,557,560,563,566,569,572,575 -o mean >> "$clean_intsv_name.$clean_bed_name.regional_bedtoolsMapMean.$cohort.tsv"
			echo "[$(date +"%Y-%m-%d %H:%M:%S")] output: $clean_intsv_name.$clean_bed_name.regional_bedtoolsMapMean.$cohort.tsv"
		fi

		if [ "$cohort" == "HBCC" ]; then
			echo "[$(date +"%Y-%m-%d %H:%M:%S")] bedtools map HBCC"
			bedtools map -a ${clean_bed_name} -b ${filtfile}.sorted.tsv -c  5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95,98,101,104,107,110,113,116,119,122,125,128,131,134,137,140,143,146,149,152,155,158,161,164,167,170,173,176,179,182,185,188,191,194,197,200,203,206,209,212,215,218,221,224,227,230,233,236,239,242,245,248,251,254,257,260,263,266,269,272,275,278,281,284,287,290,293,296,299,302,305,308,311,314,317,320,323,326,329,332,335,338,341,344,347,350,353,356,359,362,365,368,371,374,377,380,383,386,389,392,395,398,401,404,407,410,413,416,419,422,425,428,431,434,437,440,443,446,449,452,455,458,461,464 -o mean >> "$clean_intsv_name.$clean_bed_name.regional_bedtoolsMapMean.$cohort.tsv"
			echo "[$(date +"%Y-%m-%d %H:%M:%S")] output: $clean_intsv_name.$clean_bed_name.regional_bedtoolsMapMean.$cohort.tsv"
		fi

	fi
else
	echo "Directory $outdirname does not exist."
	exit 1
fi

