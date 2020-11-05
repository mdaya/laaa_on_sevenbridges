#!/bin/bash

#Set parameters
assoc_file=$1
begin_hg19_pos=$2
end_hg19_pos=$3
out_file_prefix=$4

#Get closest coordinates to being and end position specified
cat /home/analyst/locus_zoom.R | R --vanilla --args \
   $assoc_file \
   $begin_hg19_pos \
   $end_hg19_pos \
   $out_file_prefix
