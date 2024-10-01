#!/bin/bash
#SBATCH -p std               # Partition du cluster à utiliser (ici, 'std')
#SBATCH -n 1                 # Nombre de tâches (ici, 1 tâche)
#SBATCH --cpus-per-task=4    # Nombre de cœurs par tâche (ici, 4 cœurs)
#SBATCH --mem=16G            # Mémoire allouée (ici, 16 Go)
#SBATCH --time=24:00:00      
#SBATCH -J nextflow_job      # Nom du job
#SBATCH -o stdout.txt        # Fichier de sortie standard
#SBATCH -e stderr.txt        # Fichier d'erreurs
#SBATCH --mail-type=ALL      # Recevoir toutes les notifications par email


module load nextflow/24.04.4  # Charge Nextflow version 24.04.4
module load java/22.0.2
java -version
nextflow -version


# Définir l'ID du gène
GENE_NAME="ENSG00000010810"

# Définir le répertoire de travail temporaire unique pour le job
export TMPDIR=$SCRATCH/$SLURM_JOB_ID

# Créer le répertoire temporaire si nécessaire
mkdir -p $TMPDIR



# Exécuter Nextflow avec l'ID du gène en paramètre, en utilisant les ressources définies
nextflow run /home/carrelbl/PROMISE-800-years-evolution/main.nf --gene_name $GENE_NAME 

