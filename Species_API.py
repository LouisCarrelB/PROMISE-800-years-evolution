import re
import sys
import csv
import requests
from Bio import AlignIO
import pandas as pd

def get_organism_info_uniprot(uniprot_id):
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}"
    response = requests.get(url)

    if response.status_code == 200:
        protein_data = response.json()
        organism_info = protein_data.get("organism", {})
        return get_organism_name(organism_info)

    return None

def get_organism_info_ensembl(ensembl_id):
    url = f"https://rest.ensembl.org/lookup/id/{ensembl_id}?content-type=application/json"
    response = requests.get(url)

    if response.status_code == 200:
        ensembl_data = response.json()
        return ensembl_data.get("species", None)

    return None

def get_organism_name(organism_info):
    names = organism_info.get("names", [])
    for name_entry in names:
        if name_entry.get("type", "") == "scientific":
            return name_entry.get("value", "")
    return None

def extract_id(sequence_id):
    if sequence_id.startswith("ENS"):
        return sequence_id, "ensembl"
    else:
        pattern = re.compile(r"_(.+)")
        match = pattern.search(sequence_id)
        if match:
            return match.group(1), "uniprot"
        return sequence_id, "uniprot"

def get_species_from_msa(msa_file, index_word):
    try:
        alignments = AlignIO.read(msa_file, "fasta")
        species_list = []
        index_list = []

        for record in alignments:
            sequence_id, id_type = extract_id(record.id)
            if id_type == "uniprot":
                organism_name = get_organism_info_uniprot(sequence_id)
            else:
                organism_name = get_organism_info_ensembl(sequence_id)

            if organism_name:
                species_list.append(organism_name)
            else:
                species_list.append("Unknown")

            index_list.append(index_word)

        return pd.DataFrame({"Species": species_list, "Index": index_list})

    except Exception as e:
        print(f"An error occurred: {e}")
        return pd.DataFrame()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <msa_file> <output_file> <index_word>")
        sys.exit(1)

    msa_file = sys.argv[1]
    output_file = sys.argv[2]
    index_word = sys.argv[3]

    species_table = get_species_from_msa(msa_file, index_word)
    species_table.to_csv(output_file, index=False)
    print(f"Species information saved to {output_file}")
