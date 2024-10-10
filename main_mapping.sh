#!/bin/bash

# Lire le fichier gene_list.txt et soumettre un job SLURM pour chaque gène
while IFS= read -r gene; do
  [ -n "$gene" ] && sbatch --export=ALL,gene=$gene main_a3m.sh
done < gene_list.txt