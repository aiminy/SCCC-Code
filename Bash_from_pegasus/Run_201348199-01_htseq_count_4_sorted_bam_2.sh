
#samtools sort -n "201348199-01.bam_sorted.bam.bam" "201348199-01.bam_sorted.bam.bam"_sorted.bam

htseq-count -f bam -r name -s yes -i gene_name "201348199-01.bam_sorted.bam.bam" genes.gtf  > "201348199-01"_raw_count_4_sorted_bam.count

