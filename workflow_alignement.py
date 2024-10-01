#!/usr/bin/env python3

import os
import sys
import subprocess
import shutil  # Importer shutil pour copier les fichiers


def main(gene_name, base_dir):
    # Définir les répertoires de base
    gene_dir = os.path.join('/scratch/carrelbl/DATA/', gene_name)
    a3m_results_dir = os.path.join(gene_dir, 'msa_results')

    # Parcourir chaque sous-dossier de msa_results
    for root, dirs, files in os.walk(a3m_results_dir):
        for dir in dirs:
            subdir_path = os.path.join(root, dir)
            info_file = os.path.join(subdir_path, 'info.txt')

            # Vérifier l'existence du fichier info.txt
            if os.path.exists(info_file):
                # Lire le fichier info.txt pour obtenir le transcript_id
                with open(info_file, 'r') as f:
                    transcript_id = f.readline().strip()

                # Chercher le fichier .a3m correspondant dans le même répertoire
                a3m_file = None
                for file in os.listdir(subdir_path):
                    if file.endswith('.a3m'):
                        a3m_file = os.path.join(subdir_path, file)
                        break

                if a3m_file is None:
                    print(f"No .a3m file found in {subdir_path}")
                    continue


                a3m_path = a3m_file


                Alignement_cmd = [
                    'python3',
                    os.path.join(base_dir, 'Alignement_ESG.py'),
                    gene_name,
                    transcript_id,
                    'no',  # Mettre 'all' ou 'no' selon la logique métier
                    'no',  # Mettre 'redistribtuin_yes' ou 'no' selon la logique métier
                    a3m_path
                ]
                # Exécuter la commande
                print(f"Running Alignement_ESG.py for {transcript_id} with {a3m_file}")
                print(Alignement_cmd)
                result = subprocess.run(Alignement_cmd, capture_output=True, text=True)

                # Vérifier s'il y a eu une erreur lors de l'exécution
                if result.returncode != 0:
                    print(result)
                    print(f"Error running Alignement_ESG.py for {transcript_id}:")
                    print(result.stderr)
                else:
                    print(f"Successfully processed {transcript_id}")

            else:
                print(f"info.txt not found in {subdir_path}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 workflow_alignement.py <gene_name> <base_dir>")
        sys.exit(1)
    gene_name = sys.argv[1]
    base_dir = sys.argv[2]
    main(gene_name, base_dir)