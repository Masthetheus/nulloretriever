from Bio import Entrez
import pandas as pd
import xml.etree.ElementTree as ET
import re
import time
import os

Entrez.email = "matheuspedroncassol@gmail.com"
# Entrez.api_key = "SUA_CHAVE"   # opcional: sobe o limite de 3 para 10 req/s

GCA_RE = re.compile(r"GCA_\d{9}\.\d+")
ATRIBUTOS_ORIGEM = [
    "strain", "isolate", "cultivar", "ecotype",
    "isolation_source", "host", "geo_loc_name", "collection_date",
    "env_broad_scale", "env_local_scale", "env_medium",
]


def carregar_gcas(caminho):
    """Extrai apenas GCAs (versionados) do arquivo, TSV ou TXT."""
    if not os.path.exists(caminho):
        raise FileNotFoundError(caminho)

    if caminho.endswith((".tsv", ".csv")):
        sep = "\t" if caminho.endswith(".tsv") else ","
        df = pd.read_csv(caminho, sep=sep, dtype=str)
        # escolhe a coluna com mais GCAs, independente do nome
        melhor = max(df.columns,
                     key=lambda c: df[c].fillna("").str.contains(GCA_RE).sum())
        gcas = df[melhor].dropna().str.extract(f"({GCA_RE.pattern})")[0].dropna()
        gcas = gcas.tolist()
    else:
        with open(caminho, encoding="utf-8") as f:
            gcas = GCA_RE.findall(f.read())

    return list(dict.fromkeys(gcas))  # dedup preservando ordem


def entrez_retry(func, tentativas=3, **kw):
    for t in range(tentativas):
        try:
            return func(**kw)
        except Exception:
            if t == tentativas - 1:
                raise
            time.sleep(1.5 * (t + 1))


def buscar_assemblies(lote):
    """Retorna {GCA: docsummary}, apenas para GCAs do lote."""
    termo = " OR ".join(f'"{g}"[Assembly Accession]' for g in lote)
    h = entrez_retry(Entrez.esearch, db="assembly", term=termo, retmax=len(lote) * 3)
    ids = Entrez.read(h)["IdList"]
    h.close()
    if not ids:
        return {}

    h = entrez_retry(Entrez.esummary, db="assembly", id=",".join(ids), report="full")
    resumos = Entrez.read(h, validate=False)
    h.close()

    docs = resumos["DocumentSummarySet"]["DocumentSummary"]
    alvo = set(lote)
    achados = {}
    for d in docs:
        syn = d.get("Synonym", {})
        gca = syn.get("Genbank") or d.get("AssemblyAccession", "")
        if gca in alvo:               # descarta tudo que não foi pedido
            achados[gca] = d
    return achados


def extrair_strain(doc):
    bio = doc.get("Biosource", {})
    strain = ""
    for item in bio.get("InfraspeciesList", []):
        if item.get("Sub_type", "").lower() in ("strain", "isolate", "cultivar"):
            strain = item.get("Sub_value", "")
            break
    return strain, bio.get("Isolate", "")


def extrair_bioproject(doc):
    """(accession, uid) do BioProject do GenBank."""
    lst = doc.get("GB_BioProjects") or doc.get("RS_BioProjects") or []
    if lst:
        return lst[0].get("BioprojectAccn", ""), lst[0].get("BioprojectId", "")
    return "", ""


def buscar_biosamples(uids):
    """{SAMN: {atributo: valor, ...}} via efetch em lote."""
    uids = [u for u in uids if u]
    out = {}
    for i in range(0, len(uids), 100):
        h = entrez_retry(Entrez.efetch, db="biosample",
                         id=",".join(uids[i:i + 100]), retmode="xml")
        root = ET.fromstring(h.read())
        h.close()
        for bs in root.iter("BioSample"):
            acc = bs.get("accession", "")
            info = {"titulo": bs.findtext(".//Description/Title") or ""}
            outros = []
            for a in bs.findall(".//Attributes/Attribute"):
                nome = a.get("harmonized_name") or a.get("attribute_name") or ""
                val = (a.text or "").strip()
                if not val or val.lower() in ("missing", "not applicable", "na", "n/a"):
                    continue
                if nome in ATRIBUTOS_ORIGEM:
                    info[nome] = val
                else:
                    outros.append(f"{nome}={val}")
            info["outros_atributos"] = "; ".join(outros)
            out[acc] = info
        time.sleep(0.2)
    return out


_cache_bp = {}

def resumo_bioproject(uid):
    if not uid:
        return ""
    if uid in _cache_bp:
        return _cache_bp[uid]
    try:
        h = entrez_retry(Entrez.efetch, db="bioproject", id=uid, retmode="xml")
        root = ET.fromstring(h.read())
        h.close()
        title = root.findtext(".//ProjectDescr/Title") or ""
        desc = root.findtext(".//ProjectDescr/Description") or ""
        txt = f"[{title}] {desc}".strip() if title or desc else ""
        txt = re.sub(r"\s+", " ", txt)
    except Exception as e:
        txt = f"Erro: {e}"
    _cache_bp[uid] = txt
    time.sleep(0.2)
    return txt


def processar(arquivo_entrada, arquivo_saida):
    gcas = carregar_gcas(arquivo_entrada)
    print(f"{len(gcas)} GCAs únicos lidos.")

    registros, nao_achados = [], []
    for i in range(0, len(gcas), 50):
        lote = gcas[i:i + 50]
        achados = buscar_assemblies(lote)
        nao_achados += [g for g in lote if g not in achados]

        biosamples = buscar_biosamples([d.get("BioSampleId", "") for d in achados.values()])

        for gca in lote:
            d = achados.get(gca)
            if not d:
                continue
            strain, isolate = extrair_strain(d)
            bp_acc, bp_uid = extrair_bioproject(d)
            bs = biosamples.get(d.get("BioSampleAccn", ""), {})

            registros.append({
                "GCA": gca,
                "GCF": d.get("Synonym", {}).get("RefSeq", ""),
                "Organism": d.get("Organism", ""),
                "Strain": strain or bs.get("strain", ""),
                "Isolate": isolate or bs.get("isolate", ""),
                "Isolation_source": bs.get("isolation_source", ""),
                "Host": bs.get("host", ""),
                "Geo_loc": bs.get("geo_loc_name", ""),
                "Collection_date": bs.get("collection_date", ""),
                "BioSample": d.get("BioSampleAccn", ""),
                "BioSample_title": bs.get("titulo", ""),
                "BioSample_outros": bs.get("outros_atributos", ""),
                "BioProject": bp_acc,
                "Project_summary": resumo_bioproject(bp_uid),
                "Submitter": d.get("SubmitterOrganization", ""),
            })
        print(f"  lote {i // 50 + 1}: {len(achados)}/{len(lote)} encontrados")

    pd.DataFrame(registros).to_csv(arquivo_saida, index=False, encoding="utf-8-sig")
    print(f"\n{len(registros)} registros salvos em '{arquivo_saida}'.")
    if nao_achados:
        print(f"{len(nao_achados)} GCAs sem correspondência: {nao_achados[:10]}...")


if __name__ == "__main__":
    processar("scerevisiae_datasets.tsv", "metadados_leveduras_gca.csv")
