#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J ANNOVAR
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

#wc -l /scratch/projects/bbc/Project/zLi/Mutation/S2vsS6/results/mutect_output_passed.vcf.recode.vcf
#wc -l "/scratch/projects/bbc/Project/zLi/Mutation/S2vsS6/results/mutect_output_passed.vcf.recode.vcf" 
#mutect_output_passed.vcf.recode.vcf

#/scratch/projects/bbc/NGS-tools/annovar/convert2annovar.pl -format vcf4old "/scratch/projects/bbc/Project/zLi/Mutation/S2vsS6/results/mutect_output_passed.vcf.recode.vcf" > /scratch/projects/bbc/Project/zLi/Mutation/S2vsS6/results/mutect_output_passed.avinput

/scratch/projects/bbc/NGS-tools/annovar/annotate_variation.pl -out /scratch/projects/bbc/Project/zLi/Mutation/S2vsS6/results/myanno -build mm10 /scratch/projects/bbc/Project/zLi/Mutation/S2vsS6/results/mutect_output_passed.avinput /scratch/projects/bbc/NGS-tools/annovar/mm10db/

#-out /scratch/projects/bbc/Project/zLi/Mutation/S2vsS6/results/myanno -remove -protocol refGene -operation g -nastring . -vcfinput

