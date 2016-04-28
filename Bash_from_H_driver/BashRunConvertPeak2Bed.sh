#!/bin/bash   

#sed 's/"//g' temp_file_C.txt > temp_file_C_2.txt                                                                                                                                                
#sh code/BashRunConvertPeak2Bed.sh temp_file_C_2.txt

root_dir=$(pwd)
echo "$root_dir"

while read line; do


f1=`echo "$line" | awk -F"\t" '{print $1}'`
f2=`echo "$line" | awk -F"\t" '{print $2}'`
f3=`echo "$line" | awk -F"\t" '{print $3}'`
f4=$(basename "$f3")

#echo "$f1"
#echo "$f2"
#echo "$f3"

#/Peak_chip_broad/2016-03-23-2-18_S7__bam_mm_shift_75_2_peaks.broadPeak 
f5="$root_dir"/peak_chip_broad/"$f4"

num_peak=`wc -l $f5`
echo "$num_peak"

#dir_name=$(dirname "$f1")
#file_name=$(basename "$f1")
#sample_name=`echo "$file_name" | awk -F"." '{$1="",OFS="\/";print;}'`

#echo "$dir_name"
#echo "$file_name"
#echo "$sample_name"

#myBam="$file_name" myGTF="Mus_musculus.GRCm38.83.processed_chr17.gtf" myGTFDB="Mus_musculus.GRCm38.83..processed_chr1.gtf.db" myGene="$f2" runipy H3K27me3.ipynb
#myBam="$file_name" myGTF="/media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf" myGTFDB="Mus_musculus.GRCm38.83.gtf.db" myGene="$f2" runipy H3K27me3.ipynb
#ngs.plot.r -G mm10 -R tss -C /media/DATA/metaseq/metaseq-Lluis-data/"$file_name":/media/DATA/metaseq/metaseq-Lluis-data/2016-03-20-1-1_S5_.bam -O "$f2"vsInp.tss -T "$f2" -L 5000
#cut -f1  "$file_name" > "$sample_name"_gene_name.txt
#grep Hoxa "$f1" | wc -l 
#grep Hoxb "$f1" | wc -l 
#grep Hoxc "$f1" | wc -l
#grep Hoxd "$f1" | wc -l 

#f6=$(echo $f2 | tr '"' '')

#echo "$f2"
cut -f1,2,3,4 $f5 > "$root_dir"/BedFile/"$f2"_from_broadPeak.bed

done < $1