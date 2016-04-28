module load star
#STAR_2.3.0e.Linux_x86_64/STAR --genomeDir /dir/STAR/Genome --readFilesIn /dir/*.fastq.gz --readFilesCommand zcat --outFileNamePrefix SampleA --runThreadN 20
STAR_2.3.0e.Linux_x86_64/STAR --genomeDir /dir/STAR/Genome --readFilesIn /dir/*.fastq.gz --readFilesCommand zcat --outFileNamePrefix SampleA --runThreadN 20