"""
graph: Functions to store the splice graph with conservation information.
"""

import collections
import logging

from functools import reduce
import pandas as pd 
import glob
import os
from Bio import SeqIO  
import sys 


def counts_seq_added(fasta_file):
    count = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith(">") and not line.startswith(">EN"):
                count += 1
    return count

def _gene2n_transcripts(data):
    """
    Return a dictionary from gene id to number of transcripts in the gene.
    """
    return data.groupby('GeneID').apply(
        lambda df: len(set(df.TranscriptIDCluster))).to_dict()


def _edge2stats(edge2gene_list):
    """
    Return two dictionaries.

    - edge to set of genes
    - edge to dict from gene to number of transcripts in the gene with that edge
    """
    edge2genes = {}
    edge2gene_abundance = {}  # edge to transcript abundance in each gene

    for edge, gene_list in edge2gene_list.items():
        edge2genes[edge] = set(gene_list)
        edge2gene_abundance[edge] = collections.Counter(gene_list)

    return edge2genes, edge2gene_abundance


def _transcript_weighted_conservation(gene2transcriptnumber,
                                      edge2gene_abundance):
    """"
    Return the transcript weighted conservation score of each edge as a dict.
    """
    edge2trx_cons = {}
    n_genes = len(gene2transcriptnumber)
    for (edge, gene_abundance) in edge2gene_abundance.items():
        transcript_fraction_sum = 0.0
        for (gene, transcript_abundance) in gene_abundance.items():
            transcript_fraction_sum += (transcript_abundance /
                                        gene2transcriptnumber[gene])
        edge2trx_cons[edge] = transcript_fraction_sum / n_genes

    return edge2trx_cons



def nodes_and_edges2genes_and_transcripts(  # pylint: disable=too-many-locals
        data):
    """
    Return five dictionaries:

    - nodes to genes
    - edges to genes
    - nodes to transcripts
    - edges to transcripts
    - edges to the transcript weighted conservation score
    """
    gene2transcriptnumber = _gene2n_transcripts(data)

    node2genes = collections.defaultdict(set)
    edge2gene_list = collections.defaultdict(list)
    node2transcripts = collections.defaultdict(set)
    edge2transcripts = collections.defaultdict(set)
    for gene_id, gene in data.groupby('GeneID'):  # noqa pylint: disable=too-many-nested-blocks
        for transcript_id, transcript in gene.groupby('TranscriptIDCluster'):
            subexons = transcript.sort_values('S_exon_Rank')['S_exonID']

            previous = 'start'
            node2genes[previous].update({gene_id})
            node2transcripts[previous].update({transcript_id})
            for _, subexon in enumerate(subexons):
                orthologs = subexon.split('/')
                for ortholog in orthologs:
                    node2genes[ortholog].update({gene_id})
                    edge2gene_list[(previous, ortholog)].append(gene_id)
                    node2transcripts[ortholog].update({transcript_id})
                    edge2transcripts[(previous,
                                      ortholog)].update({transcript_id})
                    previous = ortholog

            node2genes['stop'].update({gene_id})
            edge2gene_list[(previous, 'stop')].append(gene_id)
            node2transcripts['stop'].update({transcript_id})
            edge2transcripts[(previous, 'stop')].update({transcript_id})

    edge2genes, edge2gene_abundance = _edge2stats(edge2gene_list)
    edge2trx_cons = _transcript_weighted_conservation(gene2transcriptnumber,
                                                      edge2gene_abundance)

    return (node2genes, edge2genes, node2transcripts, edge2transcripts,
            edge2trx_cons)



def parcourir_repertoire_msa(repertoire):
    resultats = collections.defaultdict(set)  # Utilisez un defaultdict avec des ensembles comme valeurs
    counts = collections.defaultdict(set)
    for fichier in os.listdir(repertoire):
        if fichier.endswith(".fasta") and fichier.startswith("msa_s_exon_"):
            chemin_fichier = os.path.join(repertoire, fichier)
            identifiant_exon = fichier.split("_exon_")[1].rstrip(".fasta")  # Obtenez l'identifiant de l'exon
            sequences = set()  # Utilisez un ensemble au lieu d'une liste
            for record in SeqIO.parse(chemin_fichier, "fasta"):
                nom_transcrit = record.id.split()[0][1:]  # Obtenez l'identifiant du transcrit après ">"
                sequences.add(nom_transcrit)  # Utilisez add pour ajouter à l'ensemble

            resultats[identifiant_exon] = sequences
            counts[identifiant_exon] = counts_seq_added(chemin_fichier)
    return resultats,counts



def _get_elements(node2elements):
    """
    Return a set of elements from the node to (set of) elements dictionary.

    >>> sorted(_get_elements({'a':{1, 2}, 'b':{2, 3}, 'c':{1}, 'd':{4}}))
    [1, 2, 3, 4]
    """
    return reduce(set.union, node2elements.values())


def splice_graph_gml(  # pylint: disable=too-many-locals, too-many-arguments
        filename, node2genes, edge2genes, node2transcripts,
        edge2transcripts, edge2trx_cons, s_exon_2_char,N_PROT):
    """
    Store the splice graph in GML (Graph Modelling Language) format.

    It stores the conservation level information for nodes and edges, where
    conservation is the fraction of species that has that particular feature.
    """
    if not filename.endswith('.gml'):
        filename += '.gml'
        logging.warning(
            '.gml extension was added, the splice graph will be stored at %s',
            filename)

    all_sequences = len(_get_elements(node2genes))
    n_transcripts = len(_get_elements(node2transcripts))
    n_protein = N_PROT + n_transcripts # On compte ici le nombre de séquence total utilisées pour construire le graph (la dimension de gène est abandonnée)

    with open(filename, 'w', encoding="utf-8") as gml:
        gml.write('''
        graph [
            directed 1
            id 42
            label "splice graph of s-exons"
        ''')
        node2id = {}
        node_id = 1
    
        for node, genes in node2genes.items():
            pourcentage_of_sequences_used = 100.0 * (len(genes) / all_sequences) # On regarde ici le nombre de séquence dans le MSA par rapport à la totalitée des séquences utilisées dans tout les MSAs 
            transcripts = node2transcripts[node]
            #sequences_added = n_sequences_added[node]
            #transcript_fraction_a3m = 100.0 * ((len(transcripts)+sequences_added)/ n_protein) # Nombre de transcrits dans le noeds par rapport aux transcrits total utilisés dans thoraxe original (On a ici considéré en plus que une protéine du a3m = un transcrit)
            transcript_fraction_original_thoraxe = 100.0 * (len(transcripts) / n_transcripts) # Score original données par les transcrits des 12 espèces de thoraxe 
            genes_str = ','.join(sorted(genes))
            transcripts_str = ','.join(sorted(transcripts))
            out_str = f'''
                node [
                    id {node_id}
                    label "{node}"
                    transcript_fraction {transcript_fraction_original_thoraxe}
                    genes "{genes_str}"
                    transcripts "{transcripts_str}"'''
            if s_exon_2_char and node not in ['start', 'stop']:
                out_str += f'''
                    phylosofs "{s_exon_2_char[node]}"'''
                if node[:1] != '0':
                    out_str += f'''
                    consensus "No value"'''
            out_str += '''
                ]
            '''
            gml.write(out_str)
            node2id[node] = node_id
            node_id += 1
        for edge, genes in edge2genes.items():

            source_node = node2id.get(edge[0])
            target_node = node2id.get(edge[1])

            if source_node is not None and target_node is not None:
                #conservation = 100.0 * (len(genes) / n_genes)
                #transcript_fraction_a3m = 100.0 * (len(genes)/ n_protein)
                transcripts = edge2transcripts[edge]
                transcript_fraction = 100.0 * (len(transcripts) / n_transcripts)
                transcript_weighted_conservation = edge2trx_cons[edge]
                genes_str = ','.join(sorted(genes))
                transcripts_str = ','.join(sorted(transcripts))
                gml.write(f'''
                    edge [
                        source {source_node}
                        target {target_node}
                        transcript_fraction {transcript_fraction}
                        transcript_weighted_conservation {transcript_weighted_conservation}
                        genes "{genes_str}"
                        transcripts "{transcripts_str}"
                    ]
                ''')

        gml.write('''
        ]
        ''')

    return filename





if __name__ == "__main__":
    

    if len(sys.argv) != 2:
        print("Usage: python GML_Modification.py <gene_name>")
        sys.exit(1)         
    gene_name = sys.argv[1]

    s_exons_path = f"DATA/{gene_name}/thoraxe/s_exon_table_a3m.csv"
    msa = f"DATA/{gene_name}/New_alignement/"
    a3m = True 




    s_exons = pd.read_csv(s_exons_path)
    if os.path.isdir(msa):
        (node2genes, edge2genes, node2transcripts, edge2transcripts, edge2trx_cons) = nodes_and_edges2genes_and_transcripts(s_exons)
        if a3m :
            gene_name = gene_name
            GENE = "DATA/"+gene_name + "/"
            fichiers_a3m = glob.glob(GENE + "*.a3m")    
            if fichiers_a3m:
                    fichier_a3m = fichiers_a3m[0]
            with open(fichier_a3m, 'r') as file:
                    lines = file.readlines()
            N_PROT = sum(1 for line in lines if line.startswith('>'))
            #node2genes,n_sequences_added = parcourir_repertoire_msa(f"{msa}/")
            filename = GENE + "/new_a3m.gml"
        s_exon_2_char = {}
        splice_graph_gml(filename, node2genes, edge2genes, node2transcripts,
    edge2transcripts, edge2trx_cons, s_exon_2_char,N_PROT)
