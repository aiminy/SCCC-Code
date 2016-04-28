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

fbname=$(dirname "$f" | cut -d'/' -f8)

echo "$f"

#sample_name=`echo "$f" | awk -F"." '{print $1}'`

echo "$fbname"


cat<<EOF> ~/Script_bash/Run_"$fbname"_mutect_ANNOVAR_filtered.sh
#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J ANNOVAR
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err

#wme $f)/$(basename "$f" | cut -d. -f1)_str.avinputc -l /scratch/projects/bbc/Project/zLi/Mutation/S2vsS6/results/mutect_output_passed.vcf.recode.vcf
#wc -l "$f" 
#mutect_output_passed.vcf.recode.vcf

#cut -f1-2,4- "$f" > /scratch/projects/bbc/Project/zLi/Mutation/S7vsS9/results/mutect_output_2.txt

#cut -f1-2,4- "$f" > $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3.txt

#perl /scratch/projects/bbc/NGS-tools/Convert_mut_2_vcf.pl $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3.txt $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3.vcf


#grep -v REJECT  $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3.vcf > $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3_passed.vcf  

#$annovar/convert2annovar.pl -format vcf4old $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3_passed.vcf > $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3_mutect_passed.avinput

/nethome/axy148/annovar/convert2annovar.pl -format vcf4old $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3_passed.vcf > $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3_mutect_passed.avinput

#$annovar/annotate_variation.pl -out $(dirname $f)/myanno_mutect_passed -build mm10 $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3_mutect_passed.avinput $annovar/mm10db/

/nethome/axy148/annovar/annotate_variation.pl -out $(dirname $f)/myanno_mutect_passed -build mm10 $(dirname $f)/$(basename "$f" | cut -d. -f1)_cut_c3_mutect_passed.avinput $annovar/mm10db/


#-out $(dirname $f)/myanno -remove -protocol refGene -operation g -nastring . -vcfinput

EOF
bsub -P ANNOVAR < ~/Script_bash/Run_"$fbname"_mutect_ANNOVAR_filtered.sh
done < $1
