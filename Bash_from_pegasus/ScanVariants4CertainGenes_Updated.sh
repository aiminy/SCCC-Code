#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J ANNOVAR_bash
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

#Usage: sh Script_bash/ScanVariants4CertainGenes_Updated.sh PengGenelist.txt Passed_VEP_based_annotation.txt 

rm /scratch/projects/bbc/Project/zLi/Mutation/VEP.based.variant_function.7.txt

while read line; do

f=`echo "$line"`

ff=`echo "${f,,}"`

fff=`echo "${ff^}"`

echo "$fff"

while read line; do

f_vep=`echo "$line"`

echo "$f_vep" 

fbname=$(dirname "$f_vep" | cut -d'/' -f8)

title=`head -1 "$f_vep"` 

awk -F"\t" '$1=="Location" {OFS="\t";print "pair",$0}' "$f_vep" > /scratch/projects/bbc/Project/zLi/Mutation/VEP.based.variant_function.7.title.txt

awk -F"\t" -vgene="$fff" -vpair="$fbname" '$5==gene {OFS="\t";print pair,$0}' "$f_vep" >> /scratch/projects/bbc/Project/zLi/Mutation/VEP.based.variant_function.7.txt   

done < $2

done < $1

cat /scratch/projects/bbc/Project/zLi/Mutation/VEP.based.variant_function.7.title.txt /scratch/projects/bbc/Project/zLi/Mutation/VEP.based.variant_function.7.txt > /scratch/projects/bbc/Project/zLi/Mutation/VEP.based.variant_function.with.title.txt