from funcoes_mestrado import checar_genoma, ler_arquivo_caixa_alta
from filelock import FileLock
import os

# Check all genomes file for "completude"

global_log = "logs/genomecheck.log"
lock = FileLock(global_log + ".lock")

print("Starting genome check...")
try:
    genoma = snakemake.input.genome
    log = snakemake.log[0]
    print("genoma:", genoma)
    print("log:", log)
    os.makedirs(os.path.dirname(snakemake.output.checked), exist_ok=True)
    output = snakemake.output.checked

    ler_arquivo_caixa_alta(genoma)
    if checar_genoma(genoma):
        with open(output, 'w') as f:
            f.write(f"{genoma}: OK\n")
    else:
        with lock:
            with open(global_log, 'a') as f:
                f.write(f"Genoma incompleto: {genoma}\n")
        with open(output, 'w') as f:
            f.write(f"{genoma}: INCOMPLETO\n")
except Exception as e:
    print(f"Erro: {e}")
    with lock:
        with open(global_log, 'a') as f:
            f.write(f"Erro: {e}\n")
    raise