#!/bin/bash   
                                                                                                                                                
#sh BashRunNgsPlot.sh ../metaseq/Input4Heatmap.txt

while read line; do


f1=`echo "$line"`
#f2=`echo "$line" | awk -F"\t" '{print $2}'`

echo "$f1"
#echo "$f2"

#dir_name=$(dirname "$f1")
#file_name=$(basename "$f1")
#sample_name=`echo "$file_name" | awk -F"." '{print $1}'`

#echo "$dir_name"
#echo "$file_name"
#echo "$sample_name"

#myBam="$file_name" myGTF="Mus_musculus.GRCm38.83.processed_chr17.gtf" myGTFDB="Mus_musculus.GRCm38.83..processed_chr1.gtf.db" myGene="$f2" runipy H3K27me3.ipynb
#myBam="$file_name" myGTF="/media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf" myGTFDB="Mus_musculus.GRCm38.83.gtf.db" myGene="$f2" runipy H3K27me3.ipynb
#ngs.plot.r -G mm10 -R tss -C /media/DATA/metaseq/metaseq-Lluis-data/"$file_name":/media/DATA/metaseq/metaseq-Lluis-data/2016-03-20-1-1_S5_.bam -O "$f2"vsInp.tss -T "$f2" -L 5000
#cut -f1  "$file_name" > "$sample_name"_gene_name.txt
grep Hoxa "$f1" | wc -l 
grep Hoxb "$f1" | wc -l 
grep Hoxc "$f1" | wc -l
grep Hoxd "$f1" | wc -l 
done < $1