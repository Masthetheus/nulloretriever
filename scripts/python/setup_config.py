import yaml
import os

def get_organisms_from_input():
    organisms = []
    print("Digite os nomes dos organismos (um por linha, vazio para terminar):")
    while True:
        org = input("> ").strip()
        if not org:
            break
        organisms.append(org)
    return organisms

def get_organisms_from_genomes_dir():
    genomes_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "data", "genomes"))
    files = os.listdir(genomes_dir)
    # Retorna apenas arquivos (ignora subdiretórios)
    organisms = [f for f in files if os.path.isfile(os.path.join(genomes_dir, f))]
    return organisms

def get_organisms_from_file():
    path = input("Informe o caminho do arquivo com a lista de organismos: ").strip()
    with open(path) as f:
        organisms = [line.strip() for line in f if line.strip()]
    return organisms

def main():
    print("=== Gerador de config.yaml para o pipeline ===")
    print("Escolha o modo de entrada dos organismos:")
    print("1 - Digitar manualmente")
    print("2 - Possuo um arquivo com a lista de organismos")
    print("3 - Os genomas alvos estão no diretório de genomas (data/genomes)")
    modo = input("Opção (1, 2 ou 3): ").strip()

    if modo == "1":
        organisms = get_organisms_from_input()
    elif modo == "2":
        organisms = get_organisms_from_file()
    elif modo == "3":
        organisms = get_organisms_from_genomes_dir()
        if not organisms:
            print("Nenhum organismo encontrado no diretório de genomas.")
            return
    else:
        print("Opção inválida!")
        return

    ks = input("Digite os valores de k separados por espaço (ex: 8 10 12): ").split()
    ks = [int(k) for k in ks]

    config = {
        "organisms": organisms,
        "k": ks,
        "imp": False
    }

    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    config_path = os.path.join(project_root, "config", "config.yaml")
    os.makedirs(os.path.dirname(config_path), exist_ok=True)

    with open(config_path, "w") as f:
        yaml.dump(config, f, sort_keys=False)

    print(f"\nArquivo {config_path} gerado com sucesso!")

if __name__ == "__main__":
    main()