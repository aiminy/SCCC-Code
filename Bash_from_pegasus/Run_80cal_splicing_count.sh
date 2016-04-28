#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J splicing_cal
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

awk '{if ($6 ~ "N") print}' /scratch/projects/bbc/Nimer_Cheng/80.sam > /scratch/projects/bbc/Nimer_Cheng/80.sj.txt
awk '{print $3,$4,$6}' /scratch/projects/bbc/Nimer_Cheng/80.sj.txt > /scratch/projects/bbc/Nimer_Cheng/80_tmp2
awk '{gsub("M"," M ");print}' /scratch/projects/bbc/Nimer_Cheng/80_tmp2 > /scratch/projects/bbc/Nimer_Cheng/80_tmp3
awk '{gsub("N"," N ");print}' /scratch/projects/bbc/Nimer_Cheng/80_tmp3 > /scratch/projects/bbc/Nimer_Cheng/80_tmp4
awk '{gsub("I","I ");print}' /scratch/projects/bbc/Nimer_Cheng/80_tmp4 > /scratch/projects/bbc/Nimer_Cheng/80_tmp5
awk '{gsub("D"," D ");print}' /scratch/projects/bbc/Nimer_Cheng/80_tmp5 > /scratch/projects/bbc/Nimer_Cheng/80_tmp6

awk '{for(i=1;i<=NF;i++) if ($i ~ "I") $i="";print}' /scratch/projects/bbc/Nimer_Cheng/80_tmp6 > /scratch/projects/bbc/Nimer_Cheng/80_tmp7
awk '{m=$2;for(i=3;i<=NF;i++) {if(i%2==1) m+=$i; if($i=="N") {print $1"_"m+1-$(i-1)"_"m}}}' /scratch/projects/bbc/Nimer_Cheng/80_tmp7 > /scratch/projects/bbc/Nimer_Cheng/80_tmp8
sort /scratch/projects/bbc/Nimer_Cheng/80_tmp8 | uniq > /scratch/projects/bbc/Nimer_Cheng/80_tmp9
awk '{gsub("_"," ");print}' /scratch/projects/bbc/Nimer_Cheng/80_tmp9 > /scratch/projects/bbc/Nimer_Cheng/80_tmp10
sort -k 1,1 -k 2,2n /scratch/projects/bbc/Nimer_Cheng/80_tmp10 > /scratch/projects/bbc/Nimer_Cheng/80_tmp11
awk '{print $1"_"$2"_"$3}' /scratch/projects/bbc/Nimer_Cheng/80_tmp11 > /scratch/projects/bbc/Nimer_Cheng/80_tmp12


arr=$(cat /scratch/projects/bbc/Nimer_Cheng/80_tmp12)

for i in ${arr[@]};
do
echo $i;
awk -v i=$i '{if ($1==i) print}' /scratch/projects/bbc/Nimer_Cheng/80_tmp8 > /scratch/projects/bbc/Nimer_Cheng/80_tmpfile;
c=$(wc -l /scratch/projects/bbc/Nimer_Cheng/80_tmpfile | awk '{print $1 }');
awk -v i=$i -v c=$c 'BEGIN {print i, c}' >> /scratch/projects/bbc/Nimer_Cheng/80_splicing.count;
done

