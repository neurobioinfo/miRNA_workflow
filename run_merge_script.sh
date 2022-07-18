#!/bin/bash

#SBATCH --time=03:00:00
#SBATCH --job-name=miRNA
#SBATCH --account=rrg-grouleau-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x-%j.out

module load python/3.9
module load scipy-stack

python merge_miRNA_counts_files.py "sample_list.txt" "../analysis/maturemiRNAcounts" "../analysis/genome-miRNAcounts" "../analysis"
