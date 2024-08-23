#!/bin/bash
#SBATCH -p bigm
#SBATCH -N 1
#SBATCH -n 16

# Chargement du module nécessaire
module load colabfold/1.5.2/gcc

# Arguments passés par Nextflow ou la ligne de commande
fasta_file=$1   # Chemin du fichier FASTA en entrée
output_dir=$2   # Répertoire de sortie pour les résultats

# Vérifier si les arguments sont fournis
if [ -z "$fasta_file" ] || [ -z "$output_dir" ]; then
  echo "Usage: $0 <fasta_file> <output_dir>"
  exit 1
fi

# Assurez-vous que le répertoire de sortie existe
mkdir -p $output_dir

# Définir les variables d'environnement basées sur les arguments
export X=$(dirname "$fasta_file")
export RESULT_FOLDER=$output_dir

# Exécuter colabfold_search avec les paramètres appropriés
colabfold_search --threads=32 --db-load-mode=2 "$fasta_file" "${COLABFOLD_DB}" "$RESULT_FOLDER"