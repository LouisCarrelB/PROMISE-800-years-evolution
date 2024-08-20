#!/bin/bash
#SBATCH -p bigm
#SBATCH -N 1
#SBATCH -n 16 


module load colabfold/1.5.2/gcc



export X=/scratch/carrelbl/Alt
export RESULT_FOLDER=/scratch/carrelbl/



colabfold_search  --threads=32 --db-load-mode=2 $X.fasta  ${COLABFOLD_DB} $RESULT_FOLDER
