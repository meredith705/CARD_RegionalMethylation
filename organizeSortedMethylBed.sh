#!/bin/bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
	echo "Usage: $0 <region.bed> <cohort.tsv> <filter.py> <cohortId>"
	exit 1
fi

bedregions=$1
cohorttsv=$2
filterScript=$3
cohort=$4

# get the prefix of the input tsv
clean_intsv_name=$(basename "$cohorttsv")
clean_intsv_name="${clean_intsv_name%%.*}"
echo clean in tsv name: $clean_intsv_name

# get the prefix of the input bed file
clean_bed_name=$(basename "$bedregions" .bed)
#clean_bed_name="${clean_bed_name%%.*}"
echo clean bed name: $clean_bed_name

# check that the input tsv is sorted 
if [ ! -f "$clean_intsv_name.sorted.tsv.gz" ]; then
 	# sort the input file
 	echo "[$(date +"%Y-%m-%d %H:%M:%S")] Starting sort"
 	(head -n 1 "$cohorttsv" && tail -n +2 "$cohorttsv" | sort -k 1,1 -k2,2n ) | bgzip -@16 > $clean_intsv_name.sorted.tsv.gz
 	echo "[$(date +"%Y-%m-%d %H:%M:%S")] sorted and bgzip-ed"
else
 	echo "$clean_intsv_name.sorted.tsv.gz already exists."
fi

# check that the input bed is sorted and indexed
if [ ! -f "$clean_bed_name.sorted.tsv.gz.tbi" ]; then
	echo "[$(date +"%Y-%m-%d %H:%M:%S")] Sort, bgZip and index $clean_bed_name bed regions"
	(head -n 1 $bedregions && tail -n +2 $bedregions | sort -k 1,1 -k2,2n ) | bgzip -@8 > $clean_bed_name.sorted.tsv.gz
	tabix -p bed $clean_bed_name.sorted.tsv.gz
else
	echo "$clean_bed_name.sorted.tsv.gz.tbi already sorted and indexed "
fi

# intersect the sorted input files with each other 
if [ ! -f "$clean_intsv_name.$clean_bed_name.tsv.gz" ]; then
	echo "[$(date +"%Y-%m-%d %H:%M:%S")] intersect $clean_intsv_name.sorted.tsv.gz with $clean_bed_name.sorted.tsv.gz "
	bedtools intersect -header -wa -a $clean_intsv_name.sorted.tsv.gz -b $clean_bed_name.sorted.tsv.gz | awk '!/^#/ { print | "sort -k1,1V -k2,2n"; next } { print }' | bgzip -@16 > $clean_intsv_name.$clean_bed_name.tsv.gz
else
	echo "$clean_intsv_name.$clean_bed_name.tsv.gz already intersected "
fi

# index the intersected file, maybe sort it? 
if [ ! -f "$clean_intsv_name.$clean_bed_name.tsv.gz.tbi" ]; then
	echo "[$(date +"%Y-%m-%d %H:%M:%S")] Starting index"
	tabix -p bed $clean_intsv_name.$clean_bed_name.tsv.gz
else
	echo "$clean_intsv_name.$clean_bed_name.tsv.gz.tbi already indexed "
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
		echo "[$(date +"%Y-%m-%d %H:%M:%S")] Using file: $filtfile and sorting"
		bedtools sort -header -i ${filtfile} > ${filtfile}.sorted.tsv  
		echo "[$(date +"%Y-%m-%d %H:%M:%S")] Selecting only unique positions"
		uniq ${filtfile}.sorted.tsv > ${filtfile}.sorted.uniq.tsv
		# check which cohort to run mapping 
		if [ "$cohort" == "NABEC" ]; then
			echo "[$(date +"%Y-%m-%d %H:%M:%S")] bedtools map NABEC"

			bedtools map -a ${clean_bed_name}.sorted.tsv.gz -b ${filtfile}.sorted.uniq.tsv -c 5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95,98,101,104,107,110,113,116,119,122,125,128,131,134,137,140,143,146,149,152,155,158,161,164,167,170,173,176,179,182,185,188,191,194,197,200,203,206,209,212,215,218,221,224,227,230,233,236,239,242,245,248,251,254,257,260,263,266,269,272,275,278,281,284,287,290,293,296,299,302,305,308,311,314,317,320,323,326,329,332,335,338,341,344,347,350,353,356,359,362,365,368,371,374,377,380,383,386,389,392,395,398,401,404,407,410,413,416,419,422,425,428,431,434,437,440,443,446,449,452,455,458,461,464,467,470,473,476,479,482,485,488,491,494,497,500,503,506,509,512,515,518,521,524,527,530,533,536,539,542,545,548,551,554,557,560,563,566,569,572,575 -o mean >> "$clean_intsv_name.$clean_bed_name.regional_bedtoolsMapMean.$cohort.tsv"
			echo "[$(date +"%Y-%m-%d %H:%M:%S")] output: $clean_intsv_name.$clean_bed_name.regional_bedtoolsMapMean.$cohort.tsv"
		fi

		if [ "$cohort" == "HBCC" ]; then
			echo "[$(date +"%Y-%m-%d %H:%M:%S")] bedtools map HBCC"
			bedtools map -a ${clean_bed_name}.sorted.tsv.gz -b ${filtfile}.sorted.uniq.tsv -c  5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95,98,101,104,107,110,113,116,119,122,125,128,131,134,137,140,143,146,149,152,155,158,161,164,167,170,173,176,179,182,185,188,191,194,197,200,203,206,209,212,215,218,221,224,227,230,233,236,239,242,245,248,251,254,257,260,263,266,269,272,275,278,281,284,287,290,293,296,299,302,305,308,311,314,317,320,323,326,329,332,335,338,341,344,347,350,353,356,359,362,365,368,371,374,377,380,383,386,389,392,395,398,401,404,407,410,413,416,419,422,425,428,431,434,437,440,443,446,449,452,455,458,461,464 -o mean >> "$clean_intsv_name.$clean_bed_name.regional_bedtoolsMapMean.$cohort.tsv"
			echo "[$(date +"%Y-%m-%d %H:%M:%S")] output: $clean_intsv_name.$clean_bed_name.regional_bedtoolsMapMean.$cohort.tsv"
		fi

		if [ "$cohort" == "RUSH" ]; then
			echo "[$(date +"%Y-%m-%d %H:%M:%S")] bedtools map RUSH"
			bedtools map -a ${clean_bed_name}.sorted.tsv.gz -b ${filtfile}.sorted.uniq.tsv -c  6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111,114,117,120,123,126,129,132,135,138,141,144,147,150,153,156,159,162,165,168,171,174,177,180,183,186,189,192,195,198,201,204,207,210,213,216,219,222,225,228,231,234,237,240,243,246,249,252,255,258,261,264,267,270,273,276,279,282,285,288,291,294,297,300,303,306,309,312,315,318,321,324,327,330,333,336,339,342,345,348,351,354,357,360,363,366,369,372,375,378,381,384,387,390,393,396,399,402,405,408,411,414,417,420,423,426,429,432,435,438,441,444,447,450,453,456,459,462,465,468,471,474,477,480,483,486,489,492,495,498,501,504,507,510,513,516,519,522,525,528,531,534,537,540,543,546,549,552,555,558,561,564,567,570,573,576,579,582,585,588,591,594,597,600,603,606,609,612,615,618,621,624,627,630,633,636,639,642,645,648,651,654,657,660,663,666,669,672,675,678,681,684,687,690,693,696,699,702,705,708,711,714,717,720,723,726,729,732,735,738,741,744,747,750,753,756,759,762,765,768,771,774,777,780,783,786,789,792,795,798,801,804,807,810,813,816,819,822,825,828,831,834,837,840,843,846,849,852,855,858,861,864,867,870,873,876,879,882,885,888,891,894,897,900,903,906,909,912,915,918,921,924,927,930,933,936,939,942,945,948,951,954,957,960,963,966,969,972,975,978,981,984,987,990,993,996,999,1002,1005,1008,1011,1014,1017,1020,1023,1026,1029,1032,1035,1038,1041,1044,1047 -o mean >> "$clean_intsv_name.$clean_bed_name.regional_bedtoolsMapMean.$cohort.tsv"
			echo "[$(date +"%Y-%m-%d %H:%M:%S")] output: $clean_intsv_name.$clean_bed_name.regional_bedtoolsMapMean.$cohort.tsv"
		fi


	fi
else
	echo "Directory $outdirname does not exist."
	exit 1
fi

