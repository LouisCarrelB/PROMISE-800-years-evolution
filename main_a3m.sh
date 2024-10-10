#!/bin/bash
#SBATCH -p std
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH -J mapping_$gene
#SBATCH -o stdout_a3m_$gene.txt
#SBATCH -e stderr_a3m_$gene.txt

echo "Running script for a3m for gene $gene"
python3 workflow_a3m_transcript.py $gene /home/carrelbl/PROMISE-800-years-evolution

echo "Running script for alignment for gene $gene"
python3 workflow_alignement.py $gene /home/carrelbl/PROMISE-800-years-evolution

echo "Running redistribution for gene $gene"
python3 workflow_redistribution.py $gene /home/carrelbl/PROMISE-800-years-evolution

echo "Finished processing gene $gene"