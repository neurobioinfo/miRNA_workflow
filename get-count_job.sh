#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --job-name=get-count
#SBATCH --account=rrg-grouleau-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=50M
#SBATCH --output=%x-%j.out
#SBATCH --array=1-15

echo "Starting task $SLURM_ARRAY_TASK_ID"
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_list.txt)

bash get-count_one-sample.sh $sample
