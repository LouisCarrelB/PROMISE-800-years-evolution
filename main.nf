#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Paramètres fournis via la ligne de commande ou par défaut
params.gene_name = null  // Nom du gène
params.base_dir = '/home/carrelbl/PROMISE-800-years-evolution'  // Répertoire de base contenant tous les scripts

println "Base directory: ${params.base_dir}"
println "Gene name: ${params.gene_name}"

// Générer les chemins nécessaires
params.gene_dir = "/shared/home/carrell/${params.gene_name}"        // Répertoire du gène
params.msa_results_dir = "${params.gene_dir}/msa_results/"          // Répertoire pour stocker les résultats MSA



process RunWorkflow_Alignement {
    input:
    val gene_name 

    script:
    """
    python3 ${params.base_dir}/workflow_alignement.py ${gene_name} ${params.base_dir} > script_output.log 2> script_error.log
    """
}


// Définition du workflow
workflow {

    RunWorkflow_Alignement(params.gene_name)
 
}