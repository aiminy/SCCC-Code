#Download reference genome
#mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm10.chromInfo" > mm10.genome

#! /bin/bash                                                                                                                                                                                                
#sh ~/Script_bash/BashRunConvertBamtoBigwig.sh FileBam.txt
while read line; do

#f1=`echo "$line" | awk -F"\t" '{print $1}'`                                                                                                                                                                
#f2=`echo "$line" | awk -F"\t" '{print $2}'`                                                                                                                                                                

f=`echo "$line"`

echo "$f"
sample_name=`echo "$f" | awk -F"." '{print $1}'`

echo "$sample_name"

cat > ~/Script_bash/Run_"$sample_name"_to_bw.sh <<EOF

samtools sort "$f" "$sample_name"_bam_sorted_by_position
genomeCoverageBed -ibam "$sample_name"_bam_sorted_by_position.bam -bg -g ~/Mus_musculus/mm10.genome > "$sample_name"_bam_sorted_by_position.bdg
#remove header in mm10.genome, and change name to mm10_2.genome
bedGraphToBigWig "$sample_name"_bam_sorted_by_position.bdg ~/Mus_musculus/mm10_2.genome "$sample_name"_bam_sorted_by_position_3.bw

EOF

#bsub -P Bioinformatics4count < ~/Script_bash/Run_"$sample_name"_to_bw.sh

bsub -e %J.err -P Bioinformatics4count < ~/Script_bash/Run_"$sample_name"_to_bw.sh

done < $1

#samtools sort 201348189-01.bam 201348189-01_bam_sorted_by_position
#genomeCoverageBed -ibam 201348189-01_bam_sorted_by_position.bam -bg -g mm10.genome > 201348189-01_bam_sorted_by_position.bdg
#bedtools slop -i 201348189-01_bam_sorted_by_position.bdg -g mm10.genome -b 0 
#bedGraphToBigWig 201348189-01_bam_sorted_by_position.bdg mm10_2.genome 201348189-01_bam_sorted_by_position_2.bw
#| bedClip stdin mm10.genome 201348189-01_bam_sorted_by_position.bdg.clip