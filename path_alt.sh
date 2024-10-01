#!/bin/bash
#SBATCH -p std
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH -J generate_fasta
#SBATCH -o stdout_path.txt
#SBATCH -e stderr_path.txt


python3 get_transcripts_uniq_path.py /scratch/carrelbl/DATA/ENSG00000010810/ /scratch/carrelbl/DATA/ENSG00000010810/

echo "Job finished"