#!/bin/bash

#SBATCH --job-name=snake
#SBATCH --output=snake.txt
#SBATCH --time=2-00:00:00
#SBATCH --partition=ycga
#SBATCH --nodes=1                    # number of cores and nodes
#SBATCH --cpus-per-task=16           # number of cores
#SBATCH --mem-per-cpu=10G            # shared memory, scaling with CPU request

module load BEDTools
module load SAMtools
module load miniconda
module load BWA
module load angsd
module load BCFtools 
conda activate snakemake
snakemake --cores $SLURM_CPUS_PER_TASK --scheduler greedy
