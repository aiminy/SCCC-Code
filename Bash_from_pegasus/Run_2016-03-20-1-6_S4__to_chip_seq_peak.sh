
#!/bin/bash
#BSUB -n 32
#BSUB -q parallel
#BSUB -W 05:00
#BSUB -J macs2_chip_seq
#BSUB -P macs2_chip_seq
#BSUB -o %J.out
#BSUB -e %J.err

#samtools sort "/scratch/projects/bbc/Project/Lluis/Alignment/2016-03-20-1-6_S4__BWA/2016-03-20-1-6_S4_.bam" /scratch/projects/bbc/BAM2BW/"2016-03-20-1-6_S4_"_bam_sorted_by_position
#genomeCoverageBed -ibam /scratch/projects/bbc/BAM2BW/"2016-03-20-1-6_S4_"_bam_sorted_by_position.bam -bg -g ~/Mus_musculus/mm10.genome > /scratch/projects/bbc/BAM2BW/"2016-03-20-1-6_S4_"_bam_sorted_by_position.bdg
#remove header in mm10.genome, and change name to mm10_2.genome
#LC_COLLATE=C sort -k1,1 -k2,2n /scratch/projects/bbc/BAM2BW/"2016-03-20-1-6_S4_"_bam_sorted_by_position.bdg > /scratch/projects/bbc/BAM2BW/"2016-03-20-1-6_S4_"_bam_sorted_by_position_sorted.bdg 
#bedGraphToBigWig  /scratch/projects/bbc/BAM2BW/"2016-03-20-1-6_S4_"_bam_sorted_by_position_sorted.bdg ~/Mus_musculus/mm10_2.genome /scratch/projects/bbc/BAM2BW/"2016-03-20-1-6_S4_"_bam_sorted_by_position.bw

#python /nethome/axy148/MACS/bin/macs2 callpeak -t /scratch/projects/bbc/BAM2BW/"2016-03-20-1-6_S4_"_bam_sorted_by_position.bam -c /scratch/projects/bbc/Project/Lluis/Alignment/2016-03-20-1-1_S5__BWA/2016-03-20-1-1_S5_.bam -f BAM -g mm -n "2016-03-20-1-6_S4_"_bam_mm --outdir "/scratch/projects/bbc/Peak_chip_seq/Peak_chip_broad" -B -q 0.01

#python /nethome/axy148/MACS/bin/macs2 callpeak -t /scratch/projects/bbc/BAM2BW/"2016-03-20-1-6_S4_"_bam_sorted_by_position.bam -f BAM -g mm -n "2016-03-20-1-6_S4_"_bam_mm_shift_75_2 --nomodel --shift=75 --extsize=150 --outdir "/scratch/projects/bbc/Peak_chip_seq/Peak_chip_broad" -B -q 0.1


python /nethome/axy148/MACS/bin/macs2 callpeak -t /scratch/projects/bbc/BAM2BW/"2016-03-20-1-6_S4_"_bam_sorted_by_position.bam -f BAM --broad -g mm --broad-cutoff 0.1 -n "2016-03-20-1-6_S4_"_bam_mm_shift_75_2 --nomodel --shift=75 --extsize=150 --outdir "/scratch/projects/bbc/Peak_chip_seq/Peak_chip_broad"

