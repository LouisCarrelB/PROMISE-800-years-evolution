def detect_signature(sequence):
    # Créer un dictionnaire avec les positions comme clés et les acides aminés correspondants comme valeurs
    sequence_dict = {}
    for i, aa in enumerate(sequence):
        sequence_dict[i + 1] = aa

    # Définir les positions des acides aminés dans la signature alpha
    positions_alpha = [16, 23, 13, 14, 17, 19, 21]
    # Calculer le nombre d'acides aminés corrects dans la signature alpha
    correct_alpha_count = sum(sequence_dict[pos] == "H" if pos == 16 else
                              sequence_dict[pos] == "R" if pos == 23 else
                              sequence_dict[pos] == "M" if pos == 13 else
                              sequence_dict[pos] == "V" if pos == 14 else
                              sequence_dict[pos] == "K" if pos == 17 else
                              sequence_dict[pos] == "L" if pos == 19 else
                              sequence_dict[pos] == "P" if pos == 21 else False
                              for pos in positions_alpha)
    # Calculer le score de signature alpha
    score_alpha = correct_alpha_count / len(positions_alpha) * 100

    # Définir les positions des acides aminés dans la signature beta
    positions_beta = [16, 23, 11, 15, 18]
    # Calculer le nombre d'acides aminés corrects dans la signature beta
    correct_beta_count = sum(sequence_dict[pos] == "G" if pos == 16 else
                             sequence_dict[pos] == "T" if pos == 23 else
                             sequence_dict[pos] == "G" if pos == 11 else
                             sequence_dict[pos] == "K" if pos == 15 else
                             sequence_dict[pos] == "V" if pos == 18 else False
                             for pos in positions_beta)
    # Calculer le score de signature beta
    score_beta = correct_beta_count / len(positions_beta) * 100

    return score_alpha, score_beta

sequence = "VDLWSVGCIMGEMVCHKILFPGRD-"
print(detect_signature(sequence))
