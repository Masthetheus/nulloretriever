from funcoes_mestrado import TrieBitTeste, busca_min_mutacoes

def criar_trie_exemplo():
    # Cria uma TrieBitTeste para k=4 (l=2, m=4**2=16)
    l = 2
    m = 16
    trie = TrieBitTeste(m)
    # Exemplo: insere o nulômero para v1=[0,1] (A,T), v2=3
    trie.insert([0, 1], 3)  # AT + v2=3 é nulômero
    trie.insert([0, 1], 5)  # AT + v2=5 é nulômero
    return trie

def main():
    trie = criar_trie_exemplo()
    # ATAA (v1=[0,1], v2=3) foi inserido, então NÃO é nulômero
    print("Teste 1: K-mer presente (não é nulômero)")
    print(busca_min_mutacoes(trie, "ATAA", 1))  # Esperado: None

    # ATAC (v1=[0,1], v2=1) NÃO foi inserido, então É nulômero
    print("Teste 2: K-mer ausente (é nulômero)")
    print(busca_min_mutacoes(trie, "ATAC", 1))  # Esperado: 0

    # ATAA está presente, mas ATAC está ausente. Com 1 mutação (A->C na última base), vira nulômero.
    print("Teste 3: K-mer vira nulômero com 1 mutação")
    print(busca_min_mutacoes(trie, "ATAA", 1))  # Esperado: None (pois está presente), mas se testar "ATCA" (muda T->A), pode ser 1

    # K-mer totalmente diferente, que nunca vira nulômero com 1 mutação
    print("Teste 4: K-mer nunca vira nulômero com 1 mutação")
    print(busca_min_mutacoes(trie, "GGGG", 1))  # Esperado: 0 (pois não está presente, é nulômero)

if __name__ == "__main__":
    main()