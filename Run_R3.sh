#!/bin/bash

#SBATCH -o Results/Results_%a.Rout
#SBATCH --array=1-2400
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 48:00:00

module load R/3.6

R CMD BATCH --vanilla R3_Random.R  Results/Results_${SLURM_ARRAY_TASK_ID}.Rout
