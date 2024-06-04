import requests

def get_organism_info(uniprot_id):
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}"
    try:
        response = requests.get(url, headers={"Accept": "application/json"})
        response.raise_for_status()  # Raises a HTTPError for bad responses
        protein_data = response.json()
        organism_info = protein_data.get("organism", {})
        organism_name = get_organism_name(organism_info)
        seventh_word_lineage = get_seventh_word_lineage(organism_info)
        gene_id = get_ensembl_gene_id(protein_data)  # Passing the whole protein data
        
        return organism_name, seventh_word_lineage, protein_data, gene_id
    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")
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
        return lineage[6]  # 7th word, considering 0-based indexing
    return None

def get_ensembl_gene_id(data):
    if 'dbReferences' in data:
        for ref in data['dbReferences']:
            if ref['type'] == 'Ensembl' and 'properties' in ref and 'gene ID' in ref['properties']:
                return ref['properties']['gene ID']
    return None

uniprot_id = "K9M1U5"  # Replace with your UniProt ID
organism_name, seventh_word_lineage, protein_data, gene_id = get_organism_info(uniprot_id)

if organism_name and seventh_word_lineage:
    print(f"Organism: {organism_name}")
    print(f"Seventh Word in Lineage: {seventh_word_lineage}")
    print(f"Ensembl Gene ID: {gene_id}")
    print(f"Protein Data Sample: {protein_data.get('accession')}")
else:
    print("Information not found or an error occurred.")
