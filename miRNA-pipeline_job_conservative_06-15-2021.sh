#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --job-name=miRNA
#SBATCH --account=rrg-grouleau-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x-%j.out

bash ../pratibha-miRNApipeline/final-pipeline_new.sh ../data ../analysis readSet_minusOne.txt "paired" sample_RG_minusOne.csv "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GATCGTCGGACTGTAGAACTCTGAAC"





