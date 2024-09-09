#!/usr/bin/env nextflow

params.gene_name = 'ENSG00000107643'  // Nom du gène
params.base_dir = '/scratch/carrelbl/'  // Dossier de base
params.gene_dir = "${params.base_dir}/${params.gene_name}"  // Dossier du gène
params.results_dir = "${params.gene_dir}/results/"  // Dossier pour stocker les résultats
process GenerateFasta {
    input:
    path data_files from file("${params.data_dir}/*")
    output:
    path("${params.results_dir}") into fasta_output

    script:
    """
    mkdir -p ${params.results_dir}
    python3 Fasta_for_a3m.py ${params.gene_name} ${params.base_dir}
    """
}

process GenerateMSA { 
    input:
    path fasta_files from fasta_output.collect().flatten()
    output:
    path("${params.gene_dir}/msa_results/") into msa_output

    script:
    """
    bash GenerateMSA.sh ${fasta_files} ${params.gene_dir}/msa_results/
    """

}




workflow {
    fasta_output = GenerateFasta()
    msa_output = GenerateMSA(fasta_output)
}