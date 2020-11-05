#!/bin/bash

#Set parameters
alleles_rephased_file=$1
vit_file=$2
snp_info_file=$3
sample_id_file=$4
chr=$5
begin_hg19_pos=$6
end_hg19_pos=$7

#Get closest coordinates to being and end position specified
cat /home/analyst/get_coord.R | R --vanilla --args \
   $snp_info_file \
   tmp_${begin_hg19_pos}-${end_hg19_pos}.txt \
   $begin_hg19_pos \
   $end_hg19_pos
new_begin_pos=`cut -f1 tmp_${begin_hg19_pos}-${end_hg19_pos}.txt`
new_end_pos=`cut -f2 tmp_${begin_hg19_pos}-${end_hg19_pos}.txt`

#Extract the information needed to construct the LAAA model
out_file_prefix=`basename $vit_file | sed 's/.0.Viterbi.txt//'`
out_file_prefix=${out_file_prefix}_${new_begin_pos}_${new_end_pos}_
python /home/analyst/create_dose_frames.py \
   $alleles_rephased_file \
   $vit_file \
   $snp_info_file \
   $sample_id_file \
   $new_begin_pos $new_end_pos \
   $out_file_prefix
