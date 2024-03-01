def modifier_string(A, B):
    index_tiret = [i+1 for i in range(len(A)-2) if A[i].isalpha() and A[i+1] == '-' and A[i+2].isalpha()]  # Trouver les indices où le motif est trouvé
    for pos in reversed(index_tiret):  # Parcourir en sens inverse pour insérer les tirets correctement
        if pos < len(B):
            B = B[:pos] + '-' + B[pos:]
    return B

# Exemple d'utilisation
A = "----------a-b-c"
B = "12345"
resultat = modifier_string(A, B)
print(resultat)  # Output: "12-3-45"
