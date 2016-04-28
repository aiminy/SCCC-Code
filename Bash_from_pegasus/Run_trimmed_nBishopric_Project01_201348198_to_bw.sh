
samtools sort "trimmed_nBishopric_Project01_201348198.bam" "trimmed_nBishopric_Project01_201348198"_bam_sorted_by_position
genomeCoverageBed -ibam "trimmed_nBishopric_Project01_201348198"_bam_sorted_by_position.bam -bg -g ~/Mus_musculus/mm10.genome > "trimmed_nBishopric_Project01_201348198"_bam_sorted_by_position.bdg
#remove header in mm10.genome, and change name to mm10_2.genome
bedGraphToBigWig "trimmed_nBishopric_Project01_201348198"_bam_sorted_by_position.bdg ~/Mus_musculus/mm10_2.genome "trimmed_nBishopric_Project01_201348198"_bam_sorted_by_position_3.bw

