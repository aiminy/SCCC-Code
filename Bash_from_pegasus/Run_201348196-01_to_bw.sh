
samtools sort "201348196-01.bam" "201348196-01"_bam_sorted_by_position
genomeCoverageBed -ibam "201348196-01"_bam_sorted_by_position.bam -bg -g ~/Mus_musculus/mm10.genome > "201348196-01"_bam_sorted_by_position.bdg
#remove header in mm10.genome, and change name to mm10_2.genome
bedGraphToBigWig "201348196-01"_bam_sorted_by_position.bdg ~/Mus_musculus/mm10_2.genome "201348196-01"_bam_sorted_by_position_3.bw

