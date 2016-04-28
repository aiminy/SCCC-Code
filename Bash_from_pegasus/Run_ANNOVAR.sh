#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J ANNOVAR_bash
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err


#sh Run_ANNOVAR.sh File2BeANNO.txt  

#$annovar/annotate_variation.pl -downdb -buildver mm10 -webfrom annovar refGene $annovar/mm10db/

#$annovar/annotate_variation.pl --buildver mm10 --downdb seq /scratch/projects/bbc/NGS-tools/annovar/mm10db/mm10_seq

#$annovar/retrieve_seq_from_fasta.pl $annovar/mm10db/mm10_refGene.txt -seqdir $annovar/mm10db/mm10_seq -format refGene -outfile $annovar/mm10db/mm10_refGeneMrna.fa




while read line; do

f=`echo "$line"`

fbname=$(basename "$f" | cut -d. -f1)

echo "$f"

#sample_name=`echo "$f" | awk -F"." '{print $1}'`

echo "$fbname"


cat<<EOF> ~/Script_bash/Run_"$fbname"_ANNOVAR.sh
#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J ANNOVAR
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

#wc -l /scratch/projects/bbc/Project/zLi/Mutation/S2vsS6/results/mutect_output_passed.vcf.recode.vcf
#wc -l "$f" 
#mutect_output_passed.vcf.recode.vcf

$annovar/convert2annovar.pl -format vcf4old "$f" > $(dirname $f)/$(basename "$f" | cut -d. -f1)_str.avinput

$annovar/annotate_variation.pl -out $(dirname $f)/myanno_str_new -build mm10 $(dirname $f)/$(basename "$f" | cut -d. -f1)_str.avinput $annovar/mm10db/

#-out $(dirname $f)/myanno -remove -protocol refGene -operation g -nastring . -vcfinput

EOF
bsub -P ANNOVAR < ~/Script_bash/Run_"$fbname"_ANNOVAR.sh
done < $1
