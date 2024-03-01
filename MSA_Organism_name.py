import re 
import glob
from Bio import SeqIO, AlignIO
import requests
import csv 
import pandas as pd 
import sys

def get_organism_info(uniprot_id):
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}"
    response = requests.get(url)

    if response.status_code == 200:
        protein_data = response.json()
        organism_info = protein_data.get("organism", {})
        organism_name = get_organism_name(organism_info)
        class_name = get_seventh_word_lineage(organism_info)
        embranchement_name = get_5_word_lineage(organism_info)
        regne_name = get_2_word_lineage(organism_info)
        return organism_name, class_name, embranchement_name,regne_name

    return None, None, None, None


def get_organism_name(organism_info):
    names = organism_info.get("names", [])
    for name_entry in names:
        if name_entry.get("type", "") == "scientific":
            return name_entry.get("value", "")
    return None

def get_seventh_word_lineage(organism_info):
    lineage = organism_info.get("lineage", [])
    if len(lineage) >= 7:
        return lineage[6]  
    return None

def get_5_word_lineage(organism_info):
    lineage = organism_info.get("lineage", [])
    if len(lineage) >= 5:
        return lineage[4]  
    return None

def get_2_word_lineage(organism_info):
    lineage = organism_info.get("lineage", [])
    if len(lineage) >= 2:
        return lineage[1]  
    return None

def extract_uni_id(sequence_id):
    pattern = re.compile(r"_(.+)")
    match = pattern.search(sequence_id)
    if match:
        return match.group(1)
    else:
        try:
            return sequence_id
        except ValueError:
            return None

def get_species_from_msa(msa_file):
    try:
        alignments = AlignIO.read(msa_file, "fasta")
        species_set = set()

        for record in alignments:
            uniprot_id = extract_uni_id(record.id)
            organism_name = get_organism_name(uniprot_id)

            if organism_name:
                species_set.add(organism_name)

        return list(species_set)

    except Exception as e:
        print(f"Une erreur s'est produite: {e}")
        return []

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python MSA_Organism_name.py <gene_name> ")
        sys.exit(1)         
    gene_name = sys.argv[1]
    
    GENE = "../DATA/" + gene_name
    scores_ali_can = pd.read_csv(GENE + "/inter/scores_ali_can.csv")
    scores_ali_alt = pd.read_csv(GENE + "/inter/scores_ali_alt.csv")

    scores_ali_can["UniProt_ID"] = scores_ali_can["Ligne"].apply(extract_uni_id)
    scores_ali_alt["UniProt_ID"] = scores_ali_alt["Ligne"].apply(extract_uni_id)

    organism_info_can = scores_ali_can["UniProt_ID"].apply(lambda x: get_organism_info(x))
    scores_ali_can[["Organism", "Class","Embranchement","Regne"]] = pd.DataFrame(organism_info_can.tolist(), index=scores_ali_can.index)
    scores_ali_can.fillna("No value", inplace=True)  # Fill NaN values with "NA"

    organism_info_alt = scores_ali_alt["UniProt_ID"].apply(lambda x: get_organism_info(x))
    scores_ali_alt[["Organism", "Class","Embranchement", "Regne"]] = pd.DataFrame(organism_info_alt.tolist(), index=scores_ali_alt.index)
    scores_ali_alt.fillna("No value ", inplace=True)  # Fill NaN values with "NA"

    scores_ali_can["Difference_alpha_beta"] = scores_ali_can["alpha"] - scores_ali_can["beta"]
    scores_ali_alt["Difference_alpha_beta"] = scores_ali_alt["alpha"] - scores_ali_alt["beta"]

    scores_ali_can.to_csv(GENE + "/inter/scores_ali_can_updated.csv", index=False)
    scores_ali_alt.to_csv(GENE + "/inter/scores_ali_alt_updated.csv", index=False)

    all_organisms = set()
    for organism_info in [organism_info_can, organism_info_alt]:
        for info in organism_info:
            organism_name, _, _, _ = info
            if organism_name:
                all_organisms.add(organism_name)

    with open(f"{GENE}/inter/all_organisms.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Organism"])
        for organism in all_organisms:
            writer.writerow([organism])