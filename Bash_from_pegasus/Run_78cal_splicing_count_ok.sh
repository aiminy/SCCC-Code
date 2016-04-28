#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J splicing_cal
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

#awk '{if ($6 ~ "N") print}' /scratch/projects/bbc/Nimer_Cheng/78.alignments.sam > /scratch/projects/bbc/Nimer_Cheng/78.sj.txt
#awk '{print $3,$4,$6}' /scratch/projects/bbc/Nimer_Cheng/78.sj.txt > /scratch/projects/bbc/Nimer_Cheng/78_tmp2
#awk '{gsub("M"," M ");print}' /scratch/projects/bbc/Nimer_Cheng/78_tmp2 > /scratch/projects/bbc/Nimer_Cheng/78_tmp3
#awk '{gsub("N"," N ");print}' /scratch/projects/bbc/Nimer_Cheng/78_tmp3 > /scratch/projects/bbc/Nimer_Cheng/78_tmp4
#awk '{gsub("I","I ");print}' /scratch/projects/bbc/Nimer_Cheng/78_tmp4 > /scratch/projects/bbc/Nimer_Cheng/78_tmp5
#awk '{gsub("D"," D ");print}' /scratch/projects/bbc/Nimer_Cheng/78_tmp5 > /scratch/projects/bbc/Nimer_Cheng/78_tmp6

#awk '{for(i=1;i<=NF;i++) if ($i ~ "I") $i="";print}' /scratch/projects/bbc/Nimer_Cheng/78_tmp6 > /scratch/projects/bbc/Nimer_Cheng/78_tmp7
#awk '{m=$2;for(i=3;i<=NF;i++) {if(i%2==1) m+=$i; if($i=="N") {print $1"_"m+1-$(i-1)"_"m}}}' /scratch/projects/bbc/Nimer_Cheng/78_tmp7 > /scratch/projects/bbc/Nimer_Cheng/78_tmp8
#sort /scratch/projects/bbc/Nimer_Cheng/78_tmp8 | uniq > /scratch/projects/bbc/Nimer_Cheng/78_tmp9
#awk '{gsub("_"," ");print}' /scratch/projects/bbc/Nimer_Cheng/78_tmp9 > /scratch/projects/bbc/Nimer_Cheng/78_tmp10
#sort -k 1,1 -k 2,2n /scratch/projects/bbc/Nimer_Cheng/78_tmp10 > /scratch/projects/bbc/Nimer_Cheng/78_tmp11
#awk '{print $1"_"$2"_"$3}' /scratch/projects/bbc/Nimer_Cheng/78_tmp11 > /scratch/projects/bbc/Nimer_Cheng/78_tmp12

#for i in $(cat \$(dirname $f)$(basename "$f" | cut -d. -f1)_tmp12); do echo \$i; awk -v i=\$i '{if (\$1==i) print}' $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp8 > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmpfile; c=`wc -l \$(dirname $f)$(basename "$f" | cut -d. -f1)_tmpfile | awk '{print \$1}' `; awk -v i=\$i -v c=\$c 'BEGIN {print i, c}' >> $(dirname $f)$(basename "$f" | cut -d. -f1)_splicing.count; done

#arr=( $(echo -e "a b\nc\nd e") ); for i in ${arr[@]} ;

arr=($(cat /scratch/projects/bbc/Nimer_Cheng/78_tmp12))

for i in ${arr[@]};
do
echo $i;
awk -v i=$i '{if ($1==i) print}' /scratch/projects/bbc/Nimer_Cheng/78_tmp8 > /scratch/projects/bbc/Nimer_Cheng/78_tmpfile;
c=$(wc -l /scratch/projects/bbc/Nimer_Cheng/78_tmpfile | awk '{print $1 }');
awk -v i=$i -v c=$c 'BEGIN {print i, c}' >> /scratch/projects/bbc/Nimer_Cheng/78_splicing.count;
done
