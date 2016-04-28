#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J bam2sam_each
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

samtools view -h -o /scratch/projects/bbc/Nimer_Cheng/79.sam "/scratch/projects/bbc/Nimer_Cheng/79.alignments.bam"

