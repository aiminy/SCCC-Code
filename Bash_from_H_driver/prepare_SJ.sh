awk '{if ($6 ~ "N") print}' accepted_hits.sam > tmp1
awk '{print $3,$4,$6}' tmp1 > tmp2
awk '{gsub("M"," M ");print}' tmp2 > tmp3
awk '{gsub("N"," N ");print}' tmp3 > tmp4
awk '{gsub("I","I ");print}' tmp4 > tmp5
awk '{gsub("D"," D ");print}' tmp5 > tmp6
awk '{for(i=1;i<=NF;i++) if ($i ~ "I") $i="";print}' tmp6 > tmp7
awk '{m=$2;for(i=3;i<=NF;i++) {if(i%2==1) m+=$i; if($i=="N") {print $1"_"m+1-$(i-1)"_"m}}}' tmp7 > tmp8