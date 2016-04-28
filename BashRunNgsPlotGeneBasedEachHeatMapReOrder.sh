#!/bin/bash
                                                                                                                                                
#sh BashRunNgsPlotGeneBasedEachHeatMapReOrder.sh ./File_order_list.txt

while read line; do

echo "$line"

f1=`echo "$line" | awk -F"_" '{print $4}'`
f2=`echo "$f1" | awk -F"." '{print $1}'`

echo "$f2"

#echo "$f2" > Config_Heatmap_"$f1".txt

ngs.plot.r -G mm10 -R tss -C "$line" -O heatmap_46_56_col_rm_new_"$f2"_5kb_4_26_2016 -L 5000 -RR 30 -CD 1 -CO "snow":"snow":"red2" -GO max  


done < $1



