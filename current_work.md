# REVISÃO SCRIPTS C
## 18-05 - Revisão código C.

[] main.c - kmer_extraction.c - kmerio.c
    [] main.c
        [] Rever comments
        [] Ajustar lógica pra nova abordagem
    [] kmer_extraction.c
        [] Lógica
        [] Como trata file
    [] kmerio.c
       [] Rever necessidade de manter
       [] Rever byte printing logic
[] Futuro:
    [ ] Escrever nullomer_extraction inteira em c.

## 19-05 - Revisão main.c

- Debugging fluxo:
    1. Editar direto main
    2. Traçar lógica
    3. Teste controlado
    4. Teste real
    5. Outputs

## 20-05 - Revisão main.c

- Fazer:
    [x] Revisitar Snakemake, ajustar o comando de compilação

# REVISÃO SCRIPTS PY

## 22-05 - Benchmark de extraction e statistics

[x] Gerar dados de runtime para 3 organismos diferentes, k = 14
    - GCA_000412225_2
        - extraction
              real	0m32.365s
              user	0m35.652s
              sys	0m0.992s
        - statistics
              real	1m18.303s
              user	1m17.972s
              sys	0m0.211s
    - GCA_024072835_1
        - extraction
            real	0m46.235s
            user	0m56.937s
            sys	0m1.808s
        - statistics
            real	1m11.425s
            user	1m11.123s
            sys	0m0.179s
    - GCA_033182465_1
        - extraction
            real	0m51.591s
            user	1m3.382s
            sys	0m2.099s
        - statistics
            real	1m9.927s
            user	1m9.565s
            sys	0m0.194s          
[ ] Analisar -> foco tempo




# REVISÃO SNAKEMAKE
