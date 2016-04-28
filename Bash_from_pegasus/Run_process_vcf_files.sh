#!/bin/bash

#Usage: sh ~/Script_bash/Run_process_vcf_files.sh VCF_File2BeProcessed.txt  

while read line; do

f=`echo "$line"`

fbname=$(dirname "$f" | cut -d'/' -f8)

echo "$f"


echo "$fbname"

tail -n +52 "$f" | cut -f1,2,4,5,10-11 > /scratch/projects/bbc/Project/zLi/Mutation/"$fbname"/results/Str_processed.vcf

wc -l /scratch/projects/bbc/Project/zLi/Mutation/"$fbname"/results/Str_processed.vcf

R --slave --args /scratch/projects/bbc/Project/zLi/Mutation/"$fbname"/results/Str_processed.vcf "$fbname" /scratch/projects/bbc/Project/zLi/Mutation/"$fbname"/results/Str_processed_vcf_freq_3.txt < ~/Script_bash/Run_process_one_vcf_files.R


done < $1
