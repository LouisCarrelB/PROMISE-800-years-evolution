#!/bin/bash
#SBATCH -p std
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH -J generate_fasta
#SBATCH -o stdout_a3m.txt
#SBATCH -e stderr.a3m.txt

echo "SLURM job environment initialized"
echo "Starting job"

# Activate your environment if needed
echo "Activating environment"


# Run the Python script
echo "Running script"
python3 workflow_a3m_transcript.py ENSG00000010810 /home/carrelbl/PROMISE-800-years-evolution

echo "Job finished"