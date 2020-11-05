#!/bin/bash

allele_dose_file=$1
afr_dose_file=$2
allele_afr_dose_file=$3
pheno_file=$4
r_model_file=$5
min_maf=$6
out_prefix=$7

cat /home/analyst/run_laaa.R | R --vanilla --args \
   $allele_dose_file \
   $afr_dose_file \
   $allele_afr_dose_file \
   $pheno_file \
   $r_model_file \
   $min_maf \
   $out_prefix
