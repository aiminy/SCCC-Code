#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J ANNOVAR
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

#wme /scratch/projects/bbc/Project/zLi/Mutation/S3vsS5/results/mutect_output_passed.vcf.recode.vcf)/mutect_output_passed_str.avinputc -l /scratch/projects/bbc/Project/zLi/Mutation/S2vsS6/results/mutect_output_passed.vcf.recode.vcf
wc -l "/scratch/projects/bbc/Project/zLi/Mutation/S3vsS5/results/mutect_output_passed.vcf.recode.vcf" 
#mutect_output_passed.vcf.recode.vcf

#cut -f1-2,4- "/scratch/projects/bbc/Project/zLi/Mutation/S3vsS5/results/mutect_output_passed.vcf.recode.vcf" > /scratch/projects/bbc/Project/zLi/Mutation/S7vsS9/results/mutect_output_2.txt

#cut -f1-2,4- "/scratch/projects/bbc/Project/zLi/Mutation/S3vsS5/results/mutect_output_passed.vcf.recode.vcf" > /scratch/projects/bbc/Project/zLi/Mutation/S3vsS5/results/mutect_output_passed_cut_c3.txt

#perl /scratch/projects/bbc/NGS-tools/Convert_mut_2_vcf.pl /scratch/projects/bbc/Project/zLi/Mutation/S3vsS5/results/mutect_output_passed_cut_c3.txt /scratch/projects/bbc/Project/zLi/Mutation/S3vsS5/results/mutect_output_passed_cut_c3.vcf

#/scratch/projects/bbc/NGS-tools/annovar/convert2annovar.pl -format vcf4old /scratch/projects/bbc/Project/zLi/Mutation/S3vsS5/results/mutect_output_passed_cut_c3.vcf > /scratch/projects/bbc/Project/zLi/Mutation/S3vsS5/results/mutect_output_passed_cut_c3_mutect.avinput

#/scratch/projects/bbc/NGS-tools/annovar/annotate_variation.pl -out /scratch/projects/bbc/Project/zLi/Mutation/S3vsS5/results/myanno_mutect -build mm10 /scratch/projects/bbc/Project/zLi/Mutation/S3vsS5/results/mutect_output_passed_cut_c3_mutect.avinput /scratch/projects/bbc/NGS-tools/annovar/mm10db/

#-out /scratch/projects/bbc/Project/zLi/Mutation/S3vsS5/results/myanno -remove -protocol refGene -operation g -nastring . -vcfinput

