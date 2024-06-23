import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles


def load_sequences(file_path):
    """ Charge uniquement les identifiants de séquences depuis un fichier donné.
    Chaque identifiant est situé sur une ligne commencée par '>', suivie directement par l'ID.
    """
    sequences = set()
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                # Capture l'ID jusqu'au premier espace
                identifier = line.strip().split(' ')[0][1:]  # Enlève '>' et tout après un espace
                sequences.add(identifier)
    return sequences



def analyze_changes(path_a, path_b, path_a_prime, path_b_prime):
    # Charger les séquences depuis les fichiers
    a = load_sequences(path_a)
    b = load_sequences(path_b)
    a_prime = load_sequences(path_a_prime)
    b_prime = load_sequences(path_b_prime)

    # Créer les diagrammes de Venn
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    venn2([a, a_prime], set_labels=('a', "a'"))
    venn2_circles([a, a_prime], linestyle='dashed')
    plt.title("Changements de can à can'")

    plt.subplot(1, 2, 2)
    venn2([b, b_prime], set_labels=('b', "b'"))
    venn2_circles([b, b_prime], linestyle='dashed')
    plt.title("Changements de alt à alt'")

    plt.tight_layout()
    plt.show()

    # Calculer les ajouts, pertes et échanges
    added_to_can = a_prime - a
    lost_from_can = a - a_prime
    added_to_alt = b_prime - b
    lost_from_alt = b - b_prime
    exchanged_to_can = lost_from_alt & added_to_can
    exchanged_to_alt = lost_from_can & added_to_alt

    changes_df = pd.DataFrame({
        "Added to Can'": pd.Series(list(added_to_can)),
        "Lost from Can": pd.Series(list(lost_from_can)),
        "Added to Alt'": pd.Series(list(added_to_alt)),
        "Lost from Alt": pd.Series(list(lost_from_alt)),
        "Exchanged to Can'": pd.Series(list(exchanged_to_can)),
        "Exchanged to Alt'": pd.Series(list(exchanged_to_alt))
    })

    return changes_df
