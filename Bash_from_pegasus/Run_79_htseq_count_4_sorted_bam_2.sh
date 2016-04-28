
#!/bin/bash                                                                                                                                                                                                
#BSUB -n 32                                                                                                                                                                                                
#BSUB -q general                                                                                                                                                                                           
#BSUB -W 05:00                                                                                                                                                                                             
#BSUB -J Count_BAM_each                                                                                                                                                                                    
#BSUB -P Count_BAM_each_p                                                                                                                                                                                  
#BSUB -o %J.out                                                                                                                                                                                            
#BSUB -e %J.err

#samtools sort -n "79.alignments.bam_sorted.bam.bam" "79.alignments.bam_sorted.bam.bam"_sorted.bam

htseq-count -f bam -r name -s yes -i gene_name "79.alignments.bam_sorted.bam.bam" /nethome/axy148/Mus_musculus/genes.gtf  > "79"_raw_count_4_sorted_bam.count

