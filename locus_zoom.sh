#!/bin/bash

#Set parameters
assoc_file=$1
out_file_prefix=$2
begin_hg19_pos=$3
end_hg19_pos=$4
width=$5
height=$6

#Get closest coordinates to being and end position specified
cat /home/analyst/locus_zoom.R | R --vanilla --args \
   $assoc_file \
   $out_file_prefix \
   $begin_hg19_pos \
   $end_hg19_pos \
   $width \
   $height
