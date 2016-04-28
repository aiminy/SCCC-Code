#STAR --runMode genomeGenerate --genomeDir /path/to/STAR_mm10/ --genomeFastaFiles ~/mm10_index/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --runThreadN 16
#STAR --runMode genomeGenerate --genomeDir ~/Homo_sapiens/Homo_sapiens/UCSC/hg38/Sequence/STARIndex/ --genomeFastaFiles /nethome/axy148/Homo_sapiens/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa --runThreadN 16
#/nethome/axy148/Homo_sapiens/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
STAR --genomeDir ~/Homo_sapiens/Homo_sapiens/UCSC/hg38/Sequence/STARIndex/ --readFilesIn /nethome/yxb173/Result/TrimmedData/tmp/trimmed_Koji_S1_R1.fastq /nethome/yxb173/Result/TrimmedData/tmp/trimmed_Koji_S1_R2.fastq --runThreadN 16 --outSAMtype BAM SortedByCoordinate