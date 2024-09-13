#!/usr/bin/env nextflow

// Paramètres fournis via la ligne de commande ou par défaut
params.gene_name = null  // Nom du gène
params.base_dir = '/shared/home/carrell/PROMISE-800-years-evolution'  // Répertoire de base contenant tous les scripts

println "Base directory: ${params.base_dir}"
println "Gene name: ${params.gene_name}"

// Générer les chemins nécessaires
params.gene_dir = "/shared/home/carrell/${params.gene_name}"        // Répertoire du gène
params.msa_results_dir = "${params.gene_dir}/msa_results/"          // Répertoire pour stocker les résultats MSA

process RunWorkflow_a3m_ColabFold {
    input:
    val gene_name from params.gene_name

script:
"""

python3 ${params.base_dir}/workflow_a3m_transcript.py ${gene_name} ${params.base_dir} > script_output.log 2> script_error.log
"""
}

process RunWorkflow_Alignement {
    input:
    val gene_name from params.gene_name

script:
"""
python3 ${params.base_dir}/workflow_a3m_transcript.py ${gene_name} ${params.base_dir} > script_output.log 2> script_error.log
"""
}

process Merge_MSAs {
    input:
    val gene_name from params.gene_name

script:
"""
python3 ${params.base_dir}/merge_MSA.py ${gene_name} > script_output.log 2> script_error.log
"""
}

