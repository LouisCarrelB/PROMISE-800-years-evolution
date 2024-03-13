'''
This code allows to take scores_ali_can_updated.csv as data (provided by the R script signatures.R) 
and generate data to be able to construct a phylogenetic 
tree as well as the metadata to be put in pasteML in order to annotate the same graph.
it also make a phylo tree 

'''
import pandas as pd


from Bio import Entrez, Phylo
import re 








def get_taxid(organism_name):
    Entrez.email = "louiscarrelbilliard@gmail.com" 
    handle = Entrez.esearch(db="taxonomy", term=organism_name, retmode="xml")
    record = Entrez.read(handle)
    if record["Count"] == "0":
        print(f"Aucun résultat trouvé pour {organism_name}.")
        return None
    else:
        return record["IdList"][0]


def build_tree(organism_names):
    tree = Phylo.BaseTree.Tree()
    clade = tree.clade  
    for name in organism_names:
        taxid = get_taxid(name)
        if taxid:
            child_clade = Phylo.BaseTree.Clade(name=name)  
            child_clade.taxid = taxid 
            clade.clades.append(child_clade) 
    return tree

def is_valid_organism(name):
    return bool(re.match(r'^[A-Za-z\s-]+$', name))

def make_data() : 
    # Lecture des fichiers CSV
    file_path_original = GENE + "inter/scores_ali_can_updated.csv"
    file_path_alternative = GENE + "inter/scores_ali_alt_updated.csv"

    df_original = pd.read_csv(file_path_original)
    df_alternative = pd.read_csv(file_path_alternative)
    df_original['Organism'] = df_original['Organism'].str.replace(r'\(.*?\)', '', regex=True)
    df_alternative['Organism'] = df_alternative['Organism'].str.replace(r'\(.*?\)', '', regex=True)

    # Création d'une DataFrame avec les organismes uniques des deux DataFrames
    organisms_original = set(df_original["Organism"])
    organisms_alternative = set(df_alternative["Organism"])
    all_organisms = organisms_original.union(organisms_alternative)
    # Création d'une liste pour stocker les informations sur les organismes
    organism_info = []

    # Boucle sur tous les organismes
    for organism in all_organisms:

        # Vérifie si l'organisme est présent dans les deux DataFrames, dans l'un ou dans l'autre
        if organism in organisms_original and organism in organisms_alternative:
            info = "both"
        elif organism in organisms_original:
            info = "original"
        elif organism in organisms_alternative:
            info = "alternative"
        else:
            info = None
        organism_info.append((organism, info))

    result_df = pd.DataFrame(organism_info, columns=["Organism", "Info"])
    result_df.to_csv(GENE + "inter/tempPastML_info.csv", sep="\t", index=False)


    df = pd.read_csv(GENE + "inter/tempPastML_info.csv", sep="\t")
    valid_organisms = df[df['Organism'].apply(is_valid_organism)]
    valid_organisms.to_csv(GENE + "inter/PasteML.csv", sep="\t", index=False)
    
    df = pd.read_csv(GENE + "inter/tempPastML_info.csv", sep="\t")
    organism_column = df[["Organism"]]
    #Exeption que la taxonomy de NCBI ne peux pas trouver 
    organisms_to_remove = [
        'Organism',
        'Capsaspora owczarzaki',
        'No value',
        'Cryptosporidium muris',
        'Timema poppense',
    'Pediculus humanus subsp. corporis',
    'Cryptococcus gattii serotype B',
    'Tychaedon coryphoeus'
    ]
    organism_column = organism_column[~organism_column['Organism'].isin(organisms_to_remove)]
    organism_column.to_csv(GENE + "inter/Organisms_only.csv", sep="\t", index=False)


    df = pd.read_csv(GENE + "inter/tempPastML_info.csv", sep="\t")
    filtered_df = df[~df['Organism'].isin(organisms_to_remove)]
    filtered_df.to_csv(GENE + "inter/Filtered_PasteML.csv", sep="\t", index=False)

    df = pd.read_csv(GENE + "inter/Filtered_PasteML.csv", sep="\t")
    df['Organism'] = df['Organism'].str.replace(' ', '_')
    df.to_csv(GENE + "inter/Filtered_PasteML_with_spaces.csv", sep="\t", index=False)

    return organism_column






if __name__ == "__main__": 
    gene_name = "ENSG00000107643"
    GENE = "DATA/"+ gene_name+ "/"
    organism_mat = make_data()
    organism_names = organism_mat['Organism'].tolist()
    print(organism_names)
    tree = build_tree(organism_names)
    print(tree)
    if tree.count_terminals() > 1:
        Phylo.draw(tree, do_show=True)
        Phylo.write(tree, GENE+ "inter/arbre_phylogenetique.nwk", "newick")
    else:
        print("L'arbre ne contient pas assez d'espèces pour être affiché.")

