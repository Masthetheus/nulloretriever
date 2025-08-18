from funcoes_mestrado import checar_genoma, ler_arquivo_caixa_alta
from filelock import FileLock
import os

# Check genome for ambiguous nucleotides

global_log = "logs/genomecheck.log"
lock = FileLock(global_log + ".lock")

print("Starting genome check...")
try:
    genome = snakemake.input.genome
    log = snakemake.log[0]
    print("genome:", genome)
    print("log:", log)
    os.makedirs(os.path.dirname(snakemake.output.checked), exist_ok=True)
    output = snakemake.output.checked

    ler_arquivo_caixa_alta(genome)
    if checar_genoma(genome):
        with open(output, 'w') as f:
            f.write(f"{genome}: OK\n")
    else:
        with lock:
            with open(global_log, 'a') as f:
                f.write(f"Incomplete genome: {genome}\n")
        with open(output, 'w') as f:
            f.write(f"{genome}: INCOMPLETE\n")
except Exception as e:
    print(f"Error: {e}")
    with lock:
        with open(global_log, 'a') as f:
            f.write(f"Error: {e}\n")
    raise