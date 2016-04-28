#!/bin/bash
                                                                                                                                                
#sh BashRunNgsPlotGeneBasedEachHeatMap_Cbx7_H3K27me3.sh ./Cbx7_H3K27me3.txt

mkdir Cbx7_H3K27me3_original_call_2

while read line; do

#echo "$line"

f1=`echo "$line" | awk -F"." '{print $1}'`
#f2=`echo "$f1" | awk -F"." '{print $1}'`

echo "$line"
echo "$f1"

#echo "$f2" > Config_Heatmap_"$f1".txt

ngs.plot.r -G mm10 -R tss -C "$line" -O ./Cbx7_H3K27me3_original_call_2/heatmap_max__rescale_"$f1"_5kb_4_27_2016 -L 5000 -RR 1 -CD 1 -CO "snow":"snow":"red2" -GO max

done < $1



