#!/bin/bash
#SBATCH -p cluster             # Partition du cluster
#SBATCH -n 1                   # Nombre de tâches (ici, 1 tâche)
#SBATCH --cpus-per-task=4       # Nombre de cœurs (ici, 4 cœurs)
#SBATCH --mem=16G              # Mémoire allouée (ici, 16 Go)
#SBATCH --time=01:00:00        # Temps maximum d'exécution (1 heure)
#SBATCH -J nextflow_job        # Nom du job
#SBATCH -o nextflow_output.txt # Fichier de sortie standard
#SBATCH -e nextflow_error.txt  # Fichier d'erreurs
#SBATCH --mail-type=ALL        # Notifications par email
#SBATCH --mail-user=ton.email@exemple.com  # Adresse email pour les notifications

# Charger les modules nécessaires pour Java et Nextflow
module load linux-debian10-x86_64/openjdk-11.0.8_10-gcc-8.3.0-2qkrsqu  # Charger Java
module load linux-debian10-x86_64/curl-7.76.1-gcc-8.3.0-dsz5azt  # Charger curl

# Définir l'ID du gène
GENE_NAME="ENSG00000010810"

# Exécuter Nextflow avec l'ID du gène en paramètre
NXF_VER=22.10.0 nextflow run /shared/home/carrell/PROMISE-800-years-evolution/main.nf --gene_name $GENE_NAME -resume