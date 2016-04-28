#!/bin/bash                                                                     
#BSUB -n 32                                                                     
#BSUB -q general                                                                
#BSUB -W 05:00                                                                  
#BSUB -J Sort_BAM_each                                                                 
#BSUB -P Sort_BAM_each_p                                                        
#BSUB -o %J.out                                                                 
#BSUB -e %J.err     

samtools sort -n "/scratch/projects/bbc/Nimer_Cheng/162.alignments.bam" "/scratch/projects/bbc/Nimer_Cheng/162.alignments.bam"_sorted.bam

#htseq-count -f bam -r name -s yes -i gene_name "/scratch/projects/bbc/Nimer_Cheng/162.alignments.bam"_sorted.bam genes.gtf  > "/scratch/projects/bbc/Nimer_Cheng/162"_raw_count_4_sorted_bam.txt

