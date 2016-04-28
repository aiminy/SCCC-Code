#!/bin/bash   
                                                                                                                                                
#sh BashRun_runipy_2.sh Input4Heatmap.txt

while read line; do

f1=`echo "$line" | awk -F"\t" '{print $1}'`
f2=`echo "$line" | awk -F"\t" '{print $2}'`

echo "$f1"
echo "$f2"

dir_name=$(dirname "$f1")
file_name=$(basename "$f1")
sample_name=`echo "$file_name" | awk -F"." '{print $1}'`

echo "$dir_name"
echo "$file_name"
echo "$sample_name"

#myBam="$file_name" myGTF="Mus_musculus.GRCm38.83.processed_chr17.gtf" myGTFDB="Mus_musculus.GRCm38.83..processed_chr1.gtf.db" myGene="$f2" runipy H3K27me3.ipynb
myBam="$file_name" myGTF="/media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf" myGTFDB="Mus_musculus.GRCm38.83.gtf.db" myGene="$f2" runipy H3K27me3.ipynb

done < $1


