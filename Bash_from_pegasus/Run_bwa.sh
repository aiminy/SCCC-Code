#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J myjobname
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

module load bwa
module load GATK

#java -Xmx1g -jar /share/apps/GATK/3.4.0/GenomeAnalysisTK.jar  \
#	-T RealignerTargetCreator \
#	-R /path/to/reference.fasta \
#        -o /path/to/output.intervals \
#        -known /path/to/indel_calls.vcf

#java -Xmx4g -Djava.io.tmpdir=/path/to/tmpdir \
#  -jar /share/apps/GATK/3.4.0/GenomeAnalysisTK.jar \
#  -I <lane-level.bam> \
#  -R <ref.fasta> \
#  -T IndelRealigner \
#  -targetIntervals <intervalListFromStep1Above.intervals> \
#  -o <realignedBam.bam> \
#  -known /path/to/indel_calls.vcf
#  --consensusDeterminationModel KNOWNS_ONLY \
#  -LOD 0.4

#java -jar /share/apps/picard-tools/1.103/DownsampleSam.jar

#$NGS_tools/gatk_pipeline/resources/usr/bin/vcfsorter.pl

#sh Bash_run_submit_job_3.sh FileBam.txt  
while read line; do

f=`echo "$line"`

fbname=$(basename "$f" | cut -d. -f1)

echo "$fbname"

#sample_name=`echo "$f" | awk -F"." '{print $1}'`

echo "$sample_name"

cat > ~/Script_bash/Run_"$fbname"call_var.sh <<EOF
#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J myjobname
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

java -jar /share/apps/picard-tools/1.103/DownsampleSam.jar I=$f O=$(basename "$f" | cut -d. -f1)_008.bam R=1234 P=0.08

#$NGS_tools/gatk_pipeline/resources/usr/bin/vcfsorter.pl

#java -Xmx1g -jar /share/apps/GATK/3.4.0/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R $2 \
  -o $(basename "$f" | cut -d. -f1)output.intervals \
  -known /path/to/indel_calls.vcf

java -jar /share/apps/GATK/3.4.0/GenomeAnalysisTK.jar \
   -T PrintReads \
   -R $2 \
   -I $f \
   -BQSR $(basename "$f" | cut -d. -f1).grp \
   -o $(basename "$f" | cut -d. -f1).bam

EOF
bsub -P GATK < ~/Script_bash/Run_"$fbname"call_var.sh
done < $1


