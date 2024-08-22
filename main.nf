#!/usr/bin/env nextflow

params.data_dir = 'data/'  // Chemin vers vos données d'entrée
params.result_dir = 'results/'  // Dossier pour stocker les résultats intermédiaires

process GenerateFasta {
    input:
    path data_files from file("${params.data_dir}/*") // Fichiers d'entrée
    output:
    path 'fasta_files/' into fasta_output

    script:
    """
    mkdir -p fasta_files
    python3 Fasta_for_a3m.py ${data_files}
    """
}

process GenerateMSA {
    input:
    path fasta_files from fasta_output.collect()
    output:
    path 'msa_results/' into msa_output

    script:
    """
    mkdir -p msa_results
    for fasta_file in ${fasta_files}/*.fasta; do
        bash CF_MSAgeneration.sh $fasta_file
    done
    """
}

workflow {
    fasta_output = GenerateFasta()
    msa_output = GenerateMSA(fasta_output)
    // Vous pouvez ajouter d'autres processus ici comme Alignement_ESG
}