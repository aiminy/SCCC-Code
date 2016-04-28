#!/bin/bash
                                                                                                                                                
#sh BashRunNgsPlotGeneBasedReHeatMap.sh ./27_zip_2.txt

while read line; do

#echo "$line"

f1=`echo "$line" | awk -F"\t" '{print $1}'`
f2=`echo "$line" | awk -F"\t" '{print $2}'`

#echo "$f1"

#echo "$f2" > Config_Heatmap_"$f1".txt

#ngs.plot.r -G mm10 -R tss -C Config_Heatmap_"$f1".txt -O heatmap_3_"$f1"_5kb_4_18_2016 -L 5000 -RR 1 -CD 1 -CO "blue":"snow":"red2" > output_"$f1"_heatmap_3.txt

replot.r heatmap -I "$f2" -O  re_heatmap_"$f1"_5kb_4_18_2016_max -M max

replot.r heatmap -I "$f2" -O  re_heatmap_"$f1"_5kb_4_18_2016_mean -M mean

replot.r heatmap -I "$f2" -O  re_heatmap_"$f1"_5kb_4_18_2016_window -M window

replot.r heatmap -I heatmap_2_3_5kb_4_18_2016.zip -O re_heatmap_3_5kb_4_18_2016_mean -M mean


#replot.r -I heatmap_2_1_5kb_4_18_2016.zip -O heatmap_2_1_5kb_4_18_2016_test -M mean

#dir_name=$(dirname "$f1")
#file_name=$(basename "$f1")
#sample_name=`echo "$file_name" | awk -F"." '{print $1}'`

#echo "$dir_name"
#echo "$file_name"
#echo "$sample_name"

#myBam="$file_name" myGTF="Mus_musculus.GRCm38.83.processed_chr17.gtf" myGTFDB="Mus_musculus.GRCm38.83..processed_chr1.gtf.db" myGene="$f2" runipy H3K27me3.ipynb
#myBam="$file_name" myGTF="/media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf" myGTFDB="Mus_musculus.GRCm38.83.gtf.db" myGene="$f2" runipy H3K27me3.ipynb

#
#ngs.plot.r -G mm10 -R tss -C /media/DATA/metaseq/metaseq-Lluis-data/"$file_name":/media/DATA/metaseq/metaseq-Lluis-data/2016-03-20-1-1_S5_.bam -O "$f2"vsInp.tss -T "$f2" -L 5000

#ngs.plot.r -G mm10 -R tss -C config.H3K4me3.H3K27me3.txt -O case1 -D ensembl -FL 300

#ngs.plot.r -G mm10 -R tss -C config.H3K4me3.H3K27me3.txt -O case1_5kb_rm_input -D ensembl -L 5000 -F chipseq

#ngs.plot.r -G mm10 -R tss -C config.H3K4me3.H3K27me3.txt -O case1_5kb_rm_input_2 -D ensembl -L 5000 -GO diff
#ngs.plot.r -G mm10 -R tss -C config.H3K4me3.H3K27me3.txt -O case1_5kb_rm_input_4 -L 5000 -CO "blue":"white":"red"

#ngs.plot.r -G mm10 -R tss -C config.H3K4me3.H3K27me3.txt -O case1_5kb_rm_input_5 -L 500

#ngs.plot.r -G mm10 -R tss -C config.H3K4me3.H3K27me3.txt -O case1_5kb_rm_input_6 -L 5000
#ngs.plot.r -G mm10 -R tss -C config.H3K4me3.H3K27me3.txt -O case1_5kb_rm_input_14 -L 5000 -CO "blue":"snow":"red2"

#ngs.plot.r -G mm10 -R tss -C config.H2AK119ub1.H3K27me3.txt -O case2_5kb -L 5000 -CO "blue":"snow":"red2"

#replot.r heatmap -I case2_5kb.zip -O case2_new -RR 1 -CD 1 -CO "blue":"snow":"red2" -SC global

#ngs.plot.r -G mm10 -R tss -C config.hesc.k4.txt -O hesc.k4.genebody -D ensembl -FL 300
#ngs.plot.r -G mm10 -R tss -C config.case1.txt -O case1_5kb_4_18_2016 -L 5000 -RR 1 -CD 1 -CO "blue":"snow":"red2"
#ngs.plot.r -G mm10 -R tss -C config.case2.txt -O case2_5kb_4_18_2016 -L 5000 -RR 1 -CD 1 -CO "blue":"snow":"red2"
#ngs.plot.r -G mm10 -R tss -C config.case3.txt -O case3_5kb_4_18_2016 -L 5000 -RR 1 -CD 1 -CO "blue":"snow":"red2"
#ngs.plot.r -G mm10 -R tss -C config.case4.txt -O case4_5kb_4_18_2016 -L 5000 -RR 1 -CD 1 -CO "blue":"snow":"red2"
#ngs.plot.r -G mm10 -R tss -C config.case5.txt -O case5_5kb_4_18_2016 -L 5000 -RR 1 -CD 1 -CO "blue":"snow":"red2"
#ngs.plot.r -G mm10 -R tss -C config.case6.txt -O case6_5kb_4_18_2016 -L 5000 -RR 1 -CD 1 -CO "blue":"snow":"red2"
#ngs.plot.r -G mm10 -R tss -C config.case7.txt -O case7_5kb_4_18_2016 -L 5000 -RR 1 -CD 1 -CO "blue":"snow":"red2"

#replot.r heatmap -I case2_5kb.zip -O case2_new -RR 1 -CD 1 -CO "blue":"snow":"red2" -SC global

#replot.r heatmap -I case1_5kb_4_18_2016.zip -O case1_5kb_4_18_2016_new_2 -L 5000 -RR 1 -CD 1 -CO "blue":"snow":"red2" -G0 "prod"


#ngs.plot.r -G mm10 -R tss -C config.case1-3.txt -O case1-3_5kb_4_18_2016 -L 5000 -RR 1 -CD 1 -CO "blue":"snow":"red2"
#config.case1-3.txt


done < $1



