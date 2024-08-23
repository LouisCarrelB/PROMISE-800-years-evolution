#!/bin/bash
#SBATCH -p bigm
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --job-name=nextflow_run
#SBATCH --output=nextflow_run.out

module load nextflow/20.10.0  # Chargez la version correcte de Nextflow

nextflow run main.nf --gene_name ENSG00000107643 --base_dir /scratch/carrelbl/Alt