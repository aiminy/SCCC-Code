
#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J BAM2BW
#BSUB -P BAM2BW_3_28_one
#BSUB -o %J.out
#BSUB -e %J.err

#samtools sort "/scratch/projects/bbc/Project/Lluis/Alignment/2016-03-20-1-9_S7__BWA/2016-03-20-1-9_S7_.bam" /scratch/projects/bbc/BAM2BW/"2016-03-20-1-9_S7_"_bam_sorted_by_position
#genomeCoverageBed -ibam /scratch/projects/bbc/BAM2BW/"2016-03-20-1-9_S7_"_bam_sorted_by_position.bam -bg -g ~/Mus_musculus/mm10.genome > /scratch/projects/bbc/BAM2BW/"2016-03-20-1-9_S7_"_bam_sorted_by_position.bdg
#remove header in mm10.genome, and change name to mm10_2.genome
LC_COLLATE=C sort -k1,1 -k2,2n /scratch/projects/bbc/BAM2BW/"2016-03-20-1-9_S7_"_bam_sorted_by_position.bdg > /scratch/projects/bbc/BAM2BW/"2016-03-20-1-9_S7_"_bam_sorted_by_position_sorted.bdg 
bedGraphToBigWig  /scratch/projects/bbc/BAM2BW/"2016-03-20-1-9_S7_"_bam_sorted_by_position_sorted.bdg ~/Mus_musculus/mm10_2.genome /scratch/projects/bbc/BAM2BW/"2016-03-20-1-9_S7_"_bam_sorted_by_position.bw

