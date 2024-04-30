import sys
from ete3 import Tree
import shutil

def format_name(name):
    parts = name.split('_')
    if len(parts) > 1:
        return parts[0].capitalize() + '_' + parts[1].lower()
    else:
        return name.capitalize()

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_newick_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    backup_path = file_path + ".bak"

    try:
        # Faire une copie de sauvegarde du fichier original
        shutil.copy(file_path, backup_path)
        
        # Chargement de l'arbre avec le paramètre format spécifié pour permettre des formats Newick plus complexes
        tree = Tree(open(file_path, 'r').read(), format=1)

        # Formatter chaque nom d'espèce dans l'arbre
        for leaf in tree.iter_leaves():
            leaf.name = format_name(leaf.name)

        # Réécrire l'arbre formaté dans le fichier original
        with open(file_path, 'w') as f:
            f.write(tree.write(format=1))
        
        print(f"Tree formatted and saved back to {file_path}, original saved as {backup_path}")

    except Exception as e:
        print(f"Error reading or processing the file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
