samtools view /projects/scratch/bbc/Project/Koji/Alignment/2015-11-16-21_S1_R/accepted_hits.bam | awk '{print length($10)}' | head -1000 | sort -u