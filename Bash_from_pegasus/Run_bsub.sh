#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J myjobname
#BSUB -P myprojectid
#BSUB -o %J.out
#BSUB -e %J.err
cd /nethome/axy148/Script_bash 
sh ./CheckReadLength.sh 
wait

#bsub -e %J.err -P Bioinformatics4count < $1