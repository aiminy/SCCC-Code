#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J splicing_cal_bash
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err


#sh Bash_run_submit_job_3.sh FileBam.txt  

while read line; do

f=`echo "$line"`

fbname=$(basename "$f" | cut -d. -f1)

echo "$f"

#sample_name=`echo "$f" | awk -F"." '{print $1}'`

echo "$fbname"



<<COMMENT1
####prepare the SJ positions and count of reads at SJs 

awk '{if ($6 ~ "N") print}' accepted_hits.sam > tmp1
awk '{print $3,$4,$6}' tmp1 > tmp2
awk '{gsub("M"," M ");print}' tmp2 > tmp3
awk '{gsub("N"," N ");print}' tmp3 > tmp4
awk '{gsub("I","I ");print}' tmp4 > tmp5
awk '{gsub("D"," D ");print}' tmp5 > tmp6
awk '{for(i=1;i<=NF;i++) if ($i ~ "I") $i="";print}' tmp6 > tmp7
awk '{m=$2;for(i=3;i<=NF;i++) {if(i%2==1) m+=$i; if($i=="N") {print $1"_"m+1-$(i-1)"_"m}}}' tmp7 > tmp8
######prepare the count data
sort tmp8 | uniq > tmp9
awk '{gsub("_"," ");print}' tmp9 > tmp10
sort -k 1,1 -k 2,2n tmp10 > tmp11
awk '{print $1"_"$2"_"$3}' tmp11 > tmp12
 rm splicing.count
for i in $(cat tmp12)
do
echo $i
awk -v i=$i '{if ($1==i) print}' tmp8 > tmpfile
c=`wc -l tmpfile | awk '{print $1}' `
awk -v i=$i -v c=$c 'BEGIN {print i, c}' >> splicing.count
done

#####find the corresponding gene for each splicing junctions
rm outfile1
awk '{gsub("_"," ");print}' splicing.count > tmpfile1
for i in $(cat chr)
do
awk -v i=$i '{if ($1==i) print}' tmpfile1 > tmpfile2
awk -v i=$i '{if ($1==i) print}' gene_location > tmpfile3
len=`wc -l tmpfile2 | awk '{print $1}'`
for j in `seq 1 $len`
do
m1=`awk -v j=$j '{if (NR==j) print $2}' tmpfile2`
m2=`awk -v j=$j '{if (NR==j) print $3}' tmpfile2`
awk -v m1=$m1 '{if (m1 >= $2 && m1 <= $3) print $4}' tmpfile3 > tt1
awk -v m2=$m2 '{if (m2 >= $2 && m2 <= $3) print $4}' tmpfile3 >> tt1
sort tt1 | uniq > tt2
num=`wc -l tt2  | awk '{print $1}'`
if (( $num == 0 )) 
then
awk -v j=$j '{if (NR==j) print}' tmpfile2 > tmpfile4
paste tmpfile4 tt2 >> outfile1
else
for s in `seq 1 $num`
do
awk -v j=$j '{if (NR==j) print}' tmpfile2 >> tmpfile5
done
paste tmpfile5 tt2 >> outfile1
rm tmpfile5
fi

done
done
COMMENT1


cat<<EOF> ~/Script_bash/Run_"$fbname"cal_splicing_count.sh
#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J splicing_cal
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

awk '{if (\$6 ~ "N") print}' $f > $(dirname "$f")$(basename "$f" | cut -d. -f1).sj.txt
awk '{print \$3,\$4,\$6}' $(dirname $f)$(basename "$f" | cut -d. -f1).sj.txt > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp2
awk '{gsub("M"," M ");print}' $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp2 > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp3
awk '{gsub("N"," N ");print}' $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp3 > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp4
awk '{gsub("I","I ");print}' $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp4 > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp5
awk '{gsub("D"," D ");print}' $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp5 > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp6

#awk '{for(i=1;i<=NF;i++) if ($i ~ "I") $i="";print}' tmp6 > tmp7

awk '{for(i=1;i<=NF;i++) if (\$i ~ "I") \$i="";print}' $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp6 > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp7
#awk '{m=\$2;for(i=3;i<=NF;i++) {if(i%2==1) m+=\$i; if(\$i=="N") {print \$1"_"m+1-\$(i-1)"_"m}}}' $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp7 > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp8
#sort $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp8 | uniq > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp9
#awk '{gsub("_"," ");print}' $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp9 > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp10
#sort -k 1,1 -k 2,2n $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp10 > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp11
#awk '{print \$1"_"\$2"_"\$3}' $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp11 > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp12

#for i in $(cat \$(dirname $f)$(basename "$f" | cut -d. -f1)_tmp12); do echo \$i; awk -v i=\$i '{if (\$1==i) print}' $(dirname $f)$(basename "$f" | cut -d. -f1)_tmp8 > $(dirname $f)$(basename "$f" | cut -d. -f1)_tmpfile; c=`wc -l \$(dirname $f)$(basename "$f" | cut -d. -f1)_tmpfile | awk '{print \$1}' `; awk -v i=\$i -v c=\$c 'BEGIN {print i, c}' >> $(dirname $f)$(basename "$f" | cut -d. -f1)_splicing.count; done

EOF
#bsub -P Splicing < ~/Script_bash/Run_"$fbname"cal_splicing_count.sh
done < $1
