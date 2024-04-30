from ete3 import Tree, TreeStyle, NodeStyle
import pandas as pd 
def visualize_event_from_files(newick_path, data_path):
    # Charger l'arbre à partir du fichier Newick
    with open(newick_path, 'r') as file:
        newick_str = file.read().strip()
    tree = Tree(newick_str, format=1)  # Utiliser le format NHX

    # Charger les données de l'événement
    node_data = {}
    with open(data_path, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) == 2:
                node_data[parts[0]] = bool(int(parts[1]))

    # Définir le style de l'arbre
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True

    # Définir le style de nœud pour les événements
    event_style = NodeStyle()
    event_style["bgcolor"] = "lightblue"

    no_event_style = NodeStyle()
    no_event_style["bgcolor"] = "white"

    # Appliquer le style aux nœuds selon les données
    for node in tree.traverse():
        if node.name in node_data:
            if node_data[node.name]:
                node.set_style(event_style)
            else:
                node.set_style(no_event_style)

    # Afficher l'arbre
    tree.show(tree_style=ts)

def calculate_clade_statistics(newick_path, data_path):
    # Charger l'arbre à partir du fichier Newick
    with open(newick_path, 'r') as file:
        newick_str = file.read().strip()
    tree = Tree(newick_str, format=1)  # Supposer format NHX si nécessaire

    # Charger les données de l'événement
    node_data = {}
    with open(data_path, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) == 2:
                node_data[parts[0]] = bool(int(parts[1]))

    print("Node Data:", node_data)  # Imprimer pour déboguer

    # Préparer un dictionnaire pour les statistiques des clades
    clade_stats = []

    # Calculer les statistiques pour chaque clade
    for node in tree.traverse("postorder"):  # postorder pour commencer par les feuilles
        if not node.is_leaf():
            node_event_data = [node_data.get(child.name, False) for child in node.get_children()]
            print("Node:", node.name, "Data:", node_event_data)  # Imprimer pour déboguer
            if any(node_event_data):  # Vérifier s'il y a des données à traiter
                clade_stat = {
                    'clade_root': node.name,
                    'total': len(node_event_data),
                    'event_present': sum(node_event_data),
                    'event_absent': len(node_event_data) - sum(node_event_data),
                    'percentage_event_present': 100 * sum(node_event_data) / len(node_event_data) if node_event_data else 0
                }
                clade_stats.append(clade_stat)

    # Créer un DataFrame à partir des statistiques des clades
    df = pd.DataFrame(clade_stats)
    return df





newick_path = '/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/ENSG00000107643/PASTEML_DONE/named.tree_list_species.nwk_xRSVhrp.nwk'
data_path = '/Users/louiscarrel/Documents/Alignement_Project/largescale_kinase/ENSG00000107643/PASTEML_DONE/marginal_probabilities.character_Index.model_F81.tab'


visualize_event_from_files(newick_path, data_path)


