#!/bin/bash
#SBATCH -p cluster             # Utilisation de la partition 'cluster'
#SBATCH -n 1                   # Nombre de tâches (ici, 1 tâche)
#SBATCH --cpus-per-task=4       # Nombre de cœurs (ici, 4 cœurs)
#SBATCH --mem=16G              # Mémoire allouée (ici, 16 Go)
#SBATCH --time=01:00:00        # Temps maximum d'exécution (1 heure)
#SBATCH -J nextflow_job        # Nom du job
#SBATCH -o nextflow_output.txt # Fichier de sortie standard
#SBATCH -e nextflow_error.txt  # Fichier d'erreurs
#SBATCH --mail-type=ALL        # Notifications par email pour tous les événements (début, fin, etc.)
#SBATCH --mail-user=ton.email@exemple.com  # Adresse email pour les notifications

# Charger les modules nécessaires pour Java et Nextflow
module load linux-debian10-x86_64/openjdk-11.0.8_10-gcc-8.3.0-2qkrsqu  # Charger Java
module load linux-debian10-x86_64/curl-7.76.1-gcc-8.3.0-dsz5azt

# Exécuter Nextflow directement sans srun
NXF_VER=22.10.0 nextflow run main.nf --gene_name ENSG00000107643 --base_dir ../
