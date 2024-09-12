#!/bin/bash
#SBATCH -p bigm               # Partition à utiliser (ici, bigm)
#SBATCH -n 32                 # Nombre de tâches (ici, 32 cœurs)
#SBATCH --time=1-00:00:00     # Temps maximum d'exécution (ici, 1 jour)
#SBATCH -J nextflow_job       # Nom du job
#SBATCH -o nextflow_output.txt  # Fichier de sortie standard
#SBATCH -e nextflow_error.txt   # Fichier d'erreurs
#SBATCH --mail-type=ALL       # Notifications par email pour tous les événements (début, fin, etc.)

# Charger les modules nécessaires
module load  linux-debian10-x86_64/openjdk-11.0.8_10-gcc-8.3.0-2qkrsqu
nextflow

# Exécuter le pipeline Nextflow
nextflow run main.nf --gene_name ENSG00000107643 --base_dir .