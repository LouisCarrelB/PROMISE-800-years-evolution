#!/usr/bin/env nextflow

params.gene_name = 'ENSG00000107643'  // Nom du gène
params.base_dir = '.'  // Dossier de base
params.gene_dir = "${params.base_dir}/${params.gene_name}"  // Dossier du gène
params.results_dir = "${params.gene_dir}/results_a3m/"  // Dossier pour stocker les résultats
process GenerateFasta {
    input:
    path data_files from file("${params.data_dir}/*")
    output:
    path("${params.results_dir}") into fasta_output

    script:
    ""
    mk
