#!/bin/bash

#Set parameters
assoc_file=$1
out_file_prefix=$2
chr=$3
begin_hg19_pos=$4
end_hg19_pos=$5
gene_height=$6
min_p=$7
width=$8
height=$9
base_font_size=${10}
gene_text_size=${11}

#Get closest coordinates to being and end position specified
#cat /home/analyst/locus_zoom.R | R --vanilla --args \
cat /home/analyst/locus_zoom.R | R --vanilla --args \
   $assoc_file \
   $out_file_prefix \
   $chr \
   $begin_hg19_pos \
   $end_hg19_pos \
   $gene_height \
   $min_p \
   $width \
   $height \
   $base_font_size \
   $gene_text_size 
