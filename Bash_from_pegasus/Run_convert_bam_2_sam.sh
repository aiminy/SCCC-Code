#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J bam2sam_bash
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err


#bash Run_convert_bam_2_sam.sh Cheng_bam_2.txt  

while read line; do

f=`echo "$line"`

fbname=$(basename "$f" | cut -d. -f1)

echo "$f"


echo "$fbname"

cat<<EOF> ~/Script_bash/Run_"$fbname"convert_bam2sam.sh
#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J bam2sam_each
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

samtools view -h -o $(dirname "$f")/$(basename "$f" | cut -d. -f1).sam "$f"

EOF
bsub -P Bam2sam < ~/Script_bash/Run_"$fbname"convert_bam2sam.sh
done < $1
