# Plano de Ação para Revisão de Integridade do Pipeline

**Projeto**: Nullomer Extraction Pipeline  
**Data**: 06 de Agosto de 2025  
**Objetivo**: Garantir que o pipeline atual esteja funcional e com integridade preservada

---

## 1. VERIFICAÇÃO DO AMBIENTE E DEPENDÊNCIAS

### 1.1 Ambiente Python -- OK
- **Verificar se o ambiente conda está ativo** (`nullomer-env`)
- **Validar instalação de dependências**:
  - Verificar `requirements.txt` vs pacotes instalados
  - Testar importações críticas: `numpy`, `biopython`, `bitarray`, `snakemake`, `pyyaml`
- **Executar setup se necessário**: `./setup.sh` ou instalação manual

**Comandos para validação**:
```bash
conda activate nullomer-env
pip list
python -c "import numpy, pandas, biopython, bitarray, yaml, snakemake; print('Todas as dependências OK')"
```

### 1.2 Compilação C -- OK
- **Verificar se o compilador C está disponível** (`gcc`)
- **Compilar o programa C**: `gcc scripts/c/fasta_parsing.c -o fasta_kmers_novo -O3`
- **Testar execução básica** do binário

**Comandos para validação**:
```bash
gcc --version
gcc scripts/c/fasta_parsing.c -o fasta_kmers_novo -O3
./fasta_kmers_novo --help  # ou teste básico
```

---

## 2. ESTRUTURA DE DIRETÓRIOS E CONFIGURAÇÃO

### 2.1 Verificação de Diretórios
⚠️ **ATENÇÃO CRÍTICA**: O diretório `data/` não existe atualmente

**Estrutura necessária para criar**:
```
data/
├── genomes/     # Genomas FASTA
├── raw/         # Dados brutos
└── processed/   # Dados processados
results/         # Resultados por k
logs/            # Logs de execução
benchmarks/      # Benchmarks de performance
runs/            # Status de execução
```

**Comandos para criação**:
```bash
mkdir -p data/{genomes,raw,processed}
mkdir -p results logs benchmarks runs
```

### 2.2 Configuração
- **Revisar `config/config.yaml`**:
  - Verificar se organismos listados correspondem aos arquivos FASTA disponíveis
  - Validar valores de k (atualmente apenas k=14)
  - Confirmar paths de genomas

**Pontos de verificação**:
- [ ] Organismos no config correspondem aos arquivos em `data/genomes/`
- [ ] Valores de k são adequados para os testes
- [ ] Paths estão corretos e acessíveis

---

## 3. VALIDAÇÃO DOS SCRIPTS PRINCIPAIS

### 3.1 Scripts Python Críticos
- **`funcoes_mestrado.py`** (2190 linhas): Funções core do pipeline
- **`novo_fluxo_c_e_py.py`**: Integração C+Python para cálculo de nullômeros
- **`check_genome.py`**: Validação de integridade dos genomas

### 3.2 Pontos de Atenção nos Scripts
- **Paths hardcoded**: Vários scripts contêm paths absolutos que podem não existir
- **Dependências de imports**: Verificar se todas as funções de `funcoes_mestrado.py` estão disponíveis
- **Compatibilidade de parâmetros**: Entre scripts e configuração do Snakemake

**Scripts com problemas identificados**:
- `novo_fluxo_c_e_py.py`: Path hardcoded para binário C
- `analises_finais_trie.py`: Paths absolutos para resultados
- `teste_*.py`: Múltiplos scripts com paths hardcoded

---

## 4. TESTES DE INTEGRIDADE STEP-BY-STEP

### 4.1 Teste Unitário de Componentes
1. **Testar funções core**:
   - Inicialização de trie: `inicializar_triebit_teste()`
   - Processamento de k-mers
   - Escrita de resultados: `escrever_trie_em_txt_bitarray()`

2. **Testar integração C+Python**:
   - Executar `fasta_kmers_novo` manualmente
   - Verificar output do programa C
   - Testar parsing dos k-mers em Python

**Script de teste recomendado**:
```python
# Teste básico das funções principais
from scripts.python.funcoes_mestrado import *

k = 8  # Valor pequeno para teste
l = k // 2
m = 4 ** l

# Teste de inicialização
trie = inicializar_triebit_teste(l, m)
print(f"Trie inicializada com sucesso: l={l}, m={m}")

# Teste de inserção básica
# ... adicionar testes específicos
```

### 4.2 Teste com Dataset Pequeno
- **Usar um genoma pequeno** ou subset para teste
- **Executar pipeline completo** com k pequeno (ex: k=8)
- **Verificar outputs** em cada etapa

**Configuração de teste recomendada**:
```yaml
# config/config_teste.yaml
organisms:
  - genoma_teste_pequeno
k:
  - 8
imp: false
```

---

## 5. EXECUÇÃO E VALIDAÇÃO DO SNAKEMAKE

### 5.1 Dry-run
```bash
snakemake -n -s workflow/Snakefile --configfile config/config.yaml
```

**Verificações no dry-run**:
- [ ] Todas as regras são reconhecidas
- [ ] Não há arquivos de entrada faltantes
- [ ] DAG do workflow é válido
- [ ] Paths estão corretos

### 5.2 Execução Controlada
```bash
snakemake -s workflow/Snakefile -j 1 --use-conda --configfile config/config.yaml
```

**Monitoramento durante execução**:
- [ ] Logs em tempo real sem erros críticos
- [ ] Criação de arquivos intermediários
- [ ] Uso de recursos dentro do esperado
- [ ] Progressão através das reglas do workflow

---

## 6. VERIFICAÇÕES DE SAÍDA E QUALIDADE

### 6.1 Validação de Resultados
- **Estrutura de arquivos** de saída
- **Formato dos nullômeros** (.txt files)
- **Consistência dos dados** (número esperado de nullômeros)
- **Logs de execução** sem erros críticos

**Checklist de resultados**:
- [ ] Arquivos `results/{k}/{organismo}/nullomers_{organismo}_{k}.txt` criados
- [ ] Arquivos `runs/{k}/{organismo}_done.txt` finalizados
- [ ] Logs em `logs/` sem erros críticos
- [ ] Benchmarks em `benchmarks/` com métricas válidas

### 6.2 Benchmarks e Performance
- **Uso de memória** (scripts usam `memory_profiler`)
- **Tempo de execução** por organismo/k
- **Utilização de CPU** (cores configurados)

---

## 7. PONTOS CRÍTICOS IDENTIFICADOS

### 7.1 Problemas Estruturais
- ❌ **Diretório `data/` ausente** - CRÍTICO
- ❌ **Paths absolutos hardcoded** em vários scripts - ALTO RISCO
- ❌ **Duplicação de scripts** (`workflow/scripts/` vs `scripts/`) - MÉDIO RISCO

### 7.2 Problemas de Configuração
- ⚠️ **Referência ao binário** em path incorreto no `novo_fluxo_c_e_py.py`
- ⚠️ **Dependência de paths externos** em scripts de teste
- ⚠️ **Configuração de organismos** pode não corresponder aos arquivos

### 7.3 Problemas de Integração
- ⚠️ **Interface C+Python** precisa validação
- ⚠️ **Gerenciamento de memória** para grandes genomas
- ⚠️ **Locks de arquivo** podem causar problemas em execução paralela

---

## 8. ORDEM DE EXECUÇÃO RECOMENDADA

### Fase 1: Setup Inicial
1. **Verificar e ativar ambiente Python**
2. **Compilar programa C**
3. **Criar estrutura de diretórios**

### Fase 2: Correções
4. **Corrigir paths hardcoded nos scripts**
5. **Validar configurações**
6. **Resolver duplicações de scripts**

### Fase 3: Testes
7. **Teste unitário de funções críticas**
8. **Teste com genoma pequeno**
9. **Dry-run do Snakemake**

### Fase 4: Validação
10. **Execução controlada do pipeline**
11. **Validação completa de resultados**
12. **Testes de performance**

---

## 9. CHECKLIST DE VALIDAÇÃO FINAL

### Ambiente e Dependências
- [ ] Ambiente Python (`nullomer-env`) ativo e funcional
- [ ] Todas as dependências Python instaladas
- [ ] Binário C (`fasta_kmers_novo`) compilado e testado
- [ ] Compilador GCC disponível e funcional

### Estrutura e Configuração
- [ ] Estrutura de diretórios criada
- [ ] Genomas FASTA disponíveis em `data/genomes/`
- [ ] Configuração `config.yaml` validada
- [ ] Paths corrigidos nos scripts

### Funcionalidade
- [ ] Scripts sem erros de import
- [ ] Funções principais testadas unitariamente
- [ ] Integração C+Python funcionando
- [ ] Snakemake dry-run bem-sucedido

### Execução
- [ ] Pipeline executado sem erros críticos
- [ ] Outputs gerados corretamente
- [ ] Logs sem erros críticos
- [ ] Performance dentro do esperado

### Resultados
- [ ] Arquivos de nullômeros válidos
- [ ] Estrutura de resultados correta
- [ ] Dados consistentes e esperados
- [ ] Benchmarks e métricas coletados

---

## 10. SCRIPTS DE APOIO RECOMENDADOS

### Script de Verificação Rápida
```bash
#!/bin/bash
# verificacao_rapida.sh

echo "=== Verificação Rápida do Pipeline ==="
echo "1. Verificando ambiente..."
conda list | grep -E "(numpy|pandas|biopython|snakemake)"

echo "2. Verificando estrutura..."
ls -la data/ results/ logs/ 2>/dev/null || echo "Diretórios faltando!"

echo "3. Verificando binário C..."
ls -la fasta_kmers_novo 2>/dev/null || echo "Binário C não encontrado!"

echo "4. Teste dry-run Snakemake..."
snakemake -n -s workflow/Snakefile 2>&1 | head -20
```

### Script de Setup Automático
```bash
#!/bin/bash
# setup_completo.sh

echo "=== Setup Completo do Pipeline ==="
echo "Criando estrutura de diretórios..."
mkdir -p data/{genomes,raw,processed} results logs benchmarks runs

echo "Compilando programa C..."
gcc scripts/c/fasta_parsing.c -o fasta_kmers_novo -O3

echo "Verificando dependências Python..."
pip install -r requirements.txt

echo "Setup concluído!"
```

---

## 11. CONTATOS E REFERÊNCIAS

- **Documentação principal**: `README.md`
- **Configuração**: `config/config.yaml`
- **Scripts principais**: `scripts/python/funcoes_mestrado.py`
- **Workflow**: `workflow/Snakefile`

---

**Observações Finais**:
- Este plano deve ser executado sequencialmente
- Cada falha deve ser corrigida antes de prosseguir
- Manter logs detalhados de todas as verificações
- Realizar backup antes de modificações significativas
