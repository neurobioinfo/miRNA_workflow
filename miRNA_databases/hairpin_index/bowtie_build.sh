#!/bin/bash

#SBATCH --account=rrg-grouleau-ac
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=bowtie_build
#SBATCH --output=%x-%j.out

module load StdEnv/2020 bowtie/1.3.0

bowtie-build --threads 4 -f hsa.hairpin.fa hsa.hairpin
