#mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm10.chromInfo" > mm10.genome
#samtools sort 201348189-01.bam 201348189-01_bam_sorted_by_position
#genomeCoverageBed -ibam 201348189-01_bam_sorted_by_position.bam -bg -g mm10.genome > 201348189-01_bam_sorted_by_position.bdg
#bedtools slop -i 201348189-01_bam_sorted_by_position.bdg -g mm10.genome -b 0 
bedGraphToBigWig 201348189-01_bam_sorted_by_position.bdg mm10_2.genome 201348189-01_bam_sorted_by_position_2.bw

#| bedClip stdin mm10.genome 201348189-01_bam_sorted_by_position.bdg.clip