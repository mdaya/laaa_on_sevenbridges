#!/bin/bash

#Set parameters
assoc_file=$1
out_file_prefix=$2
chr=$3
begin_hg19_pos=$4
end_hg19_pos=$5
width=$6
height=$7

#Get closest coordinates to being and end position specified
cat /home/analyst/locus_zoom.R | R --vanilla --args \
   $assoc_file \
   $out_file_prefix \
   $chr \
   $begin_hg19_pos \
   $end_hg19_pos \
   $width \
   $height
