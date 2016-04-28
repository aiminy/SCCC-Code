#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J myjobname
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

java -jar /share/apps/picard-tools/1.103/DownsampleSam.jar I=/scratch/projects/bbc/Project/zLi/Alignment/BWA/zLi_MouseExome_201519497-01_S_3_.bam O=zLi_MouseExome_201519497-01_S_3__008.bam R=1234 P=0.08

#/scratch/projects/bbc/NGS-tools/gatk_pipeline/resources/usr/bin/vcfsorter.pl

#java -Xmx1g -jar /share/apps/GATK/3.4.0/GenomeAnalysisTK.jar   -T RealignerTargetCreator   -R    -o zLi_MouseExome_201519497-01_S_3_output.intervals   -known /path/to/indel_calls.vcf

java -jar /share/apps/GATK/3.4.0/GenomeAnalysisTK.jar    -T PrintReads    -R     -I /scratch/projects/bbc/Project/zLi/Alignment/BWA/zLi_MouseExome_201519497-01_S_3_.bam    -BQSR zLi_MouseExome_201519497-01_S_3_.grp    -o zLi_MouseExome_201519497-01_S_3_.bam

