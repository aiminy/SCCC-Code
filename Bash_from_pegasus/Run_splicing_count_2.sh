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


echo "$fbname"

cat<<EOF> ~/Script_bash/Run_"$fbname"cal_splicing_count.sh
#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J splicing_cal
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

awk '{if (\$6 ~ "N") print}' $f > $(dirname "$f")/$(basename "$f" | cut -d. -f1).sj.txt
awk '{print \$3,\$4,\$6}' $(dirname $f)/$(basename "$f" | cut -d. -f1).sj.txt > $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp2
awk '{gsub("M"," M ");print}' $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp2 > $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp3
awk '{gsub("N"," N ");print}' $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp3 > $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp4
awk '{gsub("I","I ");print}' $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp4 > $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp5
awk '{gsub("D"," D ");print}' $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp5 > $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp6

awk '{for(i=1;i<=NF;i++) if (\$i ~ "I") \$i="";print}' $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp6 > $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp7
awk '{m=\$2;for(i=3;i<=NF;i++) {if(i%2==1) m+=\$i; if(\$i=="N") {print \$1"_"m+1-\$(i-1)"_"m}}}' $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp7 > $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp8
sort $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp8 | uniq > $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp9
awk '{gsub("_"," ");print}' $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp9 > $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp10
sort -k 1,1 -k 2,2n $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp10 > $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp11
awk '{print \$1"_"\$2"_"\$3}' $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp11 > $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp12


arr=\$(cat $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp12)

for i in \${arr[@]};
do
echo \$i;
awk -v i=\$i '{if (\$1==i) print}' $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmp8 > $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmpfile;
c=\$(wc -l $(dirname $f)/$(basename "$f" | cut -d. -f1)_tmpfile | awk '{print \$1 }');
awk -v i=\$i -v c=\$c 'BEGIN {print i, c}' >> $(dirname $f)/$(basename "$f" | cut -d. -f1)_splicing.count;
done

EOF
bsub -P Splicing < ~/Script_bash/Run_"$fbname"cal_splicing_count.sh
done < $1
