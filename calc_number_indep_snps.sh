#!/bin/bash

allele_dose_file=$1
laaa_results_file=$2
window_size=$3  #default 100
step_size=$4  #default 5
r2=$5  #default 0.3

#Create PLINK input file
python /home/analyst/create_tped.py $allele_dose_file

#Get list of SNPs to use
cat $laaa_results_file | cut -f1-3 | expand -t 1 | sed 's/ /:/g' > snp_list.txt

#Run PLINK
plink1.9 --tfile plink_in --extract snp_list.txt --indep-pairwise $window_size $step_size $r2
