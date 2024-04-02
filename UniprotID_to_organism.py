import requests
import pandas as pd 


def get_organism_info(uniprot_id):
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}"
    response = requests.get(url)

    if response.status_code == 200:
        protein_data = response.json()
        organism_info = protein_data.get("organism", {})
        organism_name = get_organism_name(organism_info)
        seventh_word_lineage = get_seventh_word_lineage(organism_info)
      

        return organism_name, seventh_word_lineage,protein_data

    return None, None, None
def get_organism_name(organism_info):
    names = organism_info.get("names", [])
    for name_entry in names:
        if name_entry.get("type", "") == "scientific":
            return name_entry.get("value", "")
    return None

def get_seventh_word_lineage(organism_info):
    lineage = organism_info.get("lineage", [])
    if len(lineage) >= 7:
        return lineage[6]  # 7th word, considering 0-based indexing
    return None

def gene(organism_info) :
    gene = organism_info.get("gene", [])
    return gene


uniprot_id = "A0A3Q1EHX5"  # Remplacez par votre identifiant UniProt
organism_name, seventh_word_lineage, prot = get_organism_info(uniprot_id)

if organism_name and seventh_word_lineage:
    print(f"Organism: {organism_name}")
    print(f"Seventh Word in Lineage: {seventh_word_lineage}")
else:
    print("Information not found or error occurred")

print(prot)

