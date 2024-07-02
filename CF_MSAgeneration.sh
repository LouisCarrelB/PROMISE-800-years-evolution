#!/bin/bash
#SBATCH -p bigm
#SBATCH -N 1
#SBATCH -n 16 

module load colabfold/1.5.2/gcc
colabfold_search  --threads=32 --db-load-mode=2 $X.fasta  ${COLABFOLD_DB} $RESULT_FOLDER
