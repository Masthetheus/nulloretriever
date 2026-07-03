multiplicidade_aa = {
    'A': 4, 'R': 6, 'N': 2, 'D': 2, 'C': 2,
    'E': 2, 'Q': 2, 'G': 4, 'H': 2, 'I': 3,
    'L': 6, 'K': 2, 'M': 1, 'F': 2, 'P': 4,
    'S': 6, 'T': 4, 'W': 1, 'Y': 2, 'V': 4,
    '*': 3  # Códon de parada (Stop)
}

def calcular_combinacoes(fragmento_aa):
    total_combinacoes = 1
    fragmento_aa = fragmento_aa.upper()

    for aa in fragmento_aa:
        if aa in multiplicidade_aa:
            total_combinacoes *= multiplicidade_aa[aa]
        else:
            print(f"Aviso: Aminoácido inválido encontrado -> {aa}")
            return None

    return total_combinacoes

filepath = "aas_order_all"
aas = []
with open(filepath,'r') as f:
    for line in f:
        aas.append(line.strip())

resultados = {}

for aa in aas:
    resultados[aa] = calcular_combinacoes(aa)

print(resultados)

