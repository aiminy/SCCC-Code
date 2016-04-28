samtools sort /scratch/projects/bbc/REF/GRCm38_68.fa /scratch/projects/bbc/Project/Huishi/Alignment/2015-12-09-10_S1__STAR/STAR_out.sorted.bam 2015-12-09-10_S1__STAR_STAR_out.sorted.bam 
samtools mpileup -uf /scratch/projects/bbc/REF/GRCm38_68.fa 2015-12-09-10_S1__STAR_STAR_out.sorted.bam | bcftools view -bvcg - > 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bcf
samtools mpileup -uf 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bam | bcftools view -bvcg - > 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bcf
samtools mpileup -uf 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bam | bcftools view -vcg - > 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bcf
samtools mpileup -uf 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bam | bcftools view -vcg -> 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bcf
samtools mpileup -ugf /scratch/projects/bbc/REF/GRCm38_68.fa 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bam | bcftools view -vcg - > 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.vcf
samtools mpileup -ugf /scratch/projects/bbc/REF/GRCm38_68.fa 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bam | bcftools view -vcg - > 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bcf
/share/apps/samtools/1.2/bin/samtools view 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bam | head
/share/apps/samtools/1.2/bin/samtools view 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bam | cut -f3 | head
/share/apps/samtools/1.2/bin/samtools view 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bam | cut -f3 | sort -u
/share/apps/samtools/1.2/bin/samtools view 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bam | cut -f3
/share/apps/samtools/1.2/bin/samtools view 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bam | cut -f3 | htail
/share/apps/samtools/1.2/bin/samtools view 2015-12-09-10_S1__STAR_STAR_out.sorted.bam.bam | cut -f3 | tail
 