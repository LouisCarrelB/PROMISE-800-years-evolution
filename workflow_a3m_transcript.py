#!/usr/bin/env python3

import os
import sys
import subprocess
import shutil  # Importer shutil pour copier les fichiers

def main(gene_name, base_dir):
    # Définir les répertoires de base
    gene_dir = os.path.join('/shared/home/carrell', gene_name)
    data_dir = gene_dir  # Le répertoire des données est le répertoire du gène
    results_dir1 = os.path.join(gene_dir, 'results_a3m')
    msa_results_dir = os.path.join(gene_dir, 'msa_results')

    # Assurer que les répertoires existent
    os.makedirs(results_dir1, exist_ok=True)
    os.makedirs(msa_results_dir, exist_ok=True)


    # Étape 1 : Générer les fichiers FASTA en utilisant Fasta_for_a3m.py
    print(f"Generating FASTA files for gene {gene_name}...")
    fasta_generation_cmd = [
        'python3',
        os.path.join(base_dir, 'Fasta_for_a3m.py'),
        data_dir,
        results_dir1
    ]
    result = subprocess.run(fasta_generation_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error generating FASTA files:")
        print(result.stderr)
        sys.exit(1)
    else:
        print("FASTA files generated successfully.")

    # Étape 2 : Parcourir les fichiers FASTA générés et exécuter CF_MSAgeneration.sh sur chacun
    print("Running ColabFold on generated FASTA files...")
    for root, dirs, files in os.walk(results_dir1):
        for file in files:
            if file.endswith('.fasta'):
                fasta_file = os.path.join(root, file)
                # Calculer le chemin relatif pour préserver la structure des répertoires
                rel_path = os.path.relpath(root, results_dir1)
                # Créer le répertoire de sortie correspondant
                output_subdir = os.path.join(msa_results_dir, rel_path)
                os.makedirs(output_subdir, exist_ok=True)
                # Exécuter CF_MSAgeneration.sh sur le fichier FASTA
                print(f"Processing {fasta_file}")
                colabfold_cmd = [
                    'bash',
                    os.path.join(base_dir, 'CF_MSAgeneration.sh'),
                    fasta_file,
                    output_subdir
                ]
                result = subprocess.run(colabfold_cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    print(f"Error processing {fasta_file}:")
                    print(result.stderr)
                else:
                    print(f"Successfully processed {fasta_file}")
                
                # Copie du fichier info.txt correspondant
                info_file = os.path.join(root, 'info.txt')
                if os.path.exists(info_file):
                    shutil.copy(info_file, os.path.join(output_subdir, 'info.txt'))
                    print(f"Copied info.txt to {output_subdir}")
                else:
                    print(f"No info.txt found in {root}")
            else:
                print(f"Skipping non-FASTA file: {file}")

    print("All tasks completed successfully.")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 workflow_a3m_transcript.py <gene_name> <base_dir>")
        sys.exit(1)
    gene_name = sys.argv[1]
    base_dir = sys.argv[2]
    main(gene_name, base_dir)