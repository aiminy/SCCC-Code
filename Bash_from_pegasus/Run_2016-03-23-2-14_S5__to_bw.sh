
#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J BAM2BW
#BSUB -P BAM2BW_3_28_one
#BSUB -o %J.out
#BSUB -e %J.err

samtools sort "/scratch/projects/bbc/Project/Lluis/Alignment/2016-03-23-2-14_S5__BWA/2016-03-23-2-14_S5_.bam" /scratch/projects/bbc/BAM2BW/"2016-03-23-2-14_S5_"_bam_sorted_by_position
genomeCoverageBed -ibam /scratch/projects/bbc/BAM2BW/"2016-03-23-2-14_S5_"_bam_sorted_by_position.bam -bg -g ~/Mus_musculus/mm10.genome > /scratch/projects/bbc/BAM2BW/"2016-03-23-2-14_S5_"_bam_sorted_by_position.bdg
#remove header in mm10.genome, and change name to mm10_2.genome
bedGraphToBigWig  /scratch/projects/bbc/BAM2BW/"2016-03-23-2-14_S5_"_bam_sorted_by_position.bdg ~/Mus_musculus/mm10_2.genome /scratch/projects/bbc/BAM2BW/"2016-03-23-2-14_S5_"_bam_sorted_by_position_3.bw

