#!/bin/bash
                                                                                                                                                
#sh code/BashRunProcessFiles.sh Broad_peak_18.txt Name_18_sample.txt


if [ -e "temp_file_A.txt" ]; then 
rm temp_file_A.txt;                                                                                                                                                                                     
fi

if [ -e "temp_file_B.txt" ]; then 
rm temp_file_B.txt;                                                                                                                                                                                     
fi

if [ -e "temp_file_C.txt" ]; then 
rm temp_file_C.txt;                                                                                                                                                                                     
fi

while read line; do


f1=`echo "$line"`
f2=`echo "$line" | awk -F"__" '{print $1}' | awk -F"\/" '{print $3}'`

echo "$f2\t""$f1" >> temp_file_A.txt 

#cut -f1,2,3,4 2016-03-20-1-9_S7__bam_mm_shift_75_2_peaks.broadPeak > H2AK119ub.bed

done < $1

while read line; do

f1=`echo "$line"`
f2=`echo "$line" | awk -F" " '{print $3}'`
f3=`echo "$line" | awk -F" " '{print $1}' | awk -F":" '{print $1}' | awk -F"\/" '{print $6}' | awk -F"_" '{print $1"_"$2}'`

echo "$f3\t""$f2" >>temp_file_B.txt 


#cut -f1,2,3,4 2016-03-20-1-9_S7__bam_mm_shift_75_2_peaks.broadPeak > H2AK119ub.bed

done < $2

#awk 'FILENAME==ARGV[1] {pair[$1 " " $2]; next} ($1 " " $2 in pair)' temp_file_A.txt temp_file_B.txt

#grep "$(cut -d" " -f1,2 temp_file_B.txt)" temp_file_A.txt | cut -d" " -f1,2

awk 'FNR==NR{a[$1]=$2 FS $3;next}{ print $0"\t"a[$1]}' temp_file_A.txt temp_file_B.txt | sort -uk1 > temp_file_C.txt
#awk 'FNR==NR{a[$1]=$2 FS $3;next}{ print $0}' temp_file_A.txt temp_file_B.txt | sort -uk1 > temp_file_C.txt