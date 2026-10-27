#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07 — REDE PPI CENTRADA EM NQO1 + VIAS KEGG E REACTOME
===============================================================================

O que este script faz
---------------------
1) Baixa do STRING a rede de vizinhança direta do NQO1 (Homo sapiens,
   score combinado >= 700, ou seja, alta confiança).
2) Consulta o KEGG para obter as vias que contêm NQO1 e os genes dessas vias.
3) Consulta o Reactome para obter as vias que contêm NQO1 e sua hierarquia.
4) Cruza os parceiros do STRING com as vias KEGG/Reactome.
5) Anota os nós com o log2FoldChange (DCM vs CTRL) do resultado do limma.
6) Gera figuras (rede, barras, heatmap, mapas de via) e tabelas CSV.

Saídas principais
-----------------
output/figures/nqo1_rede_ppi_string.png
output/figures/nqo1_kegg_pathways.png
output/figures/nqo1_reactome_pathways.png
output/figures/nqo1_heatmap_vias_parceiros.png
output/figures/nqo1_kegg_diagram_hsa05208.png
output/figures/nqo1_reactome_diagram_nfe2l2.png
output/tables/nqo1_string_ppi_arestas.csv
output/tables/nqo1_string_ppi_nos.csv
output/tables/nqo1_kegg_pathways.csv
output/tables/nqo1_reactome_pathways.csv
output/tables/nqo1_reactome_hierarquia.csv
output/tables/nqo1_resumo_metricas.csv

Como rodar
----------
    python tutoriais/07_PPI_NQO1_string_kegg_reactome.py

Para forçar a atualização dos dados baixados:
    NQO1_REFRESH=1 python tutoriais/07_PPI_NQO1_string_kegg_reactome.py

Observação: KEGG e Reactome são bases de vias, não de interações físicas.
A rede PPI em si vem do STRING; KEGG/Reactome são usados para anotar e
contextualizar os parceiros.
"""

from __future__ import annotations

import io
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import requests
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image, ImageDraw, ImageFont

# ---------------------------------------------------------------------------
# 0) CONFIGURAÇÃO
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
FIG = OUT / "figures"
TAB = OUT / "tables"
TMP = OUT / "temp"
CACHE = TMP / "nqo1_ppi_multidb"
for pasta in (FIG, TAB, TMP, CACHE):
    pasta.mkdir(parents=True, exist_ok=True)

GENE = "NQO1"
UNIPROT = "P15559"
SPECIES = 9606
STRING_SCORE = 700
STRING_LIMIT = 500
REFRESH = os.environ.get("NQO1_REFRESH", "0") == "1"

# Alguns símbolos do STRING têm sufixo/isoforma. Para cruzar com KEGG/Reactome,
# usamos um alias conservador quando o símbolo "limpo" existe.
ALIAS = {"SOD2-2": "SOD2"}

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "nqo1-ppi-analysis/1.0 (educational)"})


def log(msg: str) -> None:
    print(f"[NQO1] {msg}", flush=True)


# ---------------------------------------------------------------------------
# 1) FUNÇÕES DE DOWNLOAD COM CACHE
# ---------------------------------------------------------------------------
def _get(url: str, **kwargs) -> requests.Response:
    last = None
    for tentativa in range(3):
        try:
            r = SESSION.get(url, timeout=60, **kwargs)
            r.raise_for_status()
            return r
        except Exception as e:  # pragma: no cover - depende da rede
            last = e
            if tentativa < 2:
                time.sleep(1.5 * (tentativa + 1))
    raise RuntimeError(f"Falha ao acessar {url}: {last}")


def cached_text(nome: str, url: str, params: dict | None = None) -> str:
    caminho = CACHE / nome
    if caminho.exists() and not REFRESH:
        return caminho.read_text(encoding="utf-8")
    r = _get(url, params=params)
    r.encoding = "utf-8"
    caminho.write_text(r.text, encoding="utf-8")
    return r.text


def cached_json(nome: str, url: str, params: dict | None = None):
    caminho = CACHE / nome
    if caminho.exists() and not REFRESH:
        return json.loads(caminho.read_text(encoding="utf-8"))
    r = _get(url, params=params)
    dados = r.json()
    caminho.write_text(json.dumps(dados, ensure_ascii=False, indent=2), encoding="utf-8")
    return dados


def cached_bytes(nome: str, url: str) -> bytes:
    caminho = CACHE / nome
    if caminho.exists() and not REFRESH:
        return caminho.read_bytes()
    r = _get(url)
    caminho.write_bytes(r.content)
    return r.content


# ---------------------------------------------------------------------------
# 2) STRING — REDE CENTRADA EM NQO1
# ---------------------------------------------------------------------------
def baixar_rede_string() -> pd.DataFrame:
    log("Baixando rede do STRING...")
    texto = cached_text(
        "string_nqo1_network.tsv",
        "https://string-db.org/api/tsv/network",
        params={
            "identifiers": GENE,
            "species": SPECIES,
            "required_score": STRING_SCORE,
            "caller_identity": "nqo1_ppi_multidb",
            "limit": STRING_LIMIT,
        },
    )
    df = pd.read_csv(io.StringIO(texto), sep="\t")
    if df.empty:
        raise RuntimeError("O STRING não retornou interações para NQO1.")

    df = df.rename(
        columns={
            "preferredName_A": "gene_a",
            "preferredName_B": "gene_b",
            "stringId_A": "string_id_a",
            "stringId_B": "string_id_b",
        }
    )
    # Remove duplicatas A-B / B-A e guarda o maior score.
    df["par"] = df.apply(lambda r: tuple(sorted((r["gene_a"], r["gene_b"]))), axis=1)
    df = df.sort_values("score", ascending=False).drop_duplicates("par").drop(columns="par")
    df = df.reset_index(drop=True)

    df["is_nqo1_edge"] = (df["gene_a"] == GENE) | (df["gene_b"] == GENE)
    df["score_1000"] = (df["score"] * 1000).round().astype(int)

    colunas = [
        "gene_a", "gene_b", "score", "score_1000", "is_nqo1_edge",
        "nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore",
        "string_id_a", "string_id_b",
    ]
    df = df[colunas]
    log(f"STRING: {df.shape[0]} arestas únicas; "
        f"{df['is_nqo1_edge'].sum()} ligadas diretamente a NQO1.")
    return df


def montar_grafo(arestas: pd.DataFrame) -> nx.Graph:
    g = nx.Graph()
    for _, r in arestas.iterrows():
        g.add_edge(
            r["gene_a"],
            r["gene_b"],
            score=float(r["score"]),
            score_1000=int(r["score_1000"]),
            is_nqo1_edge=bool(r["is_nqo1_edge"]),
        )
    return g


# ---------------------------------------------------------------------------
# 3) EXPRESSÃO DIFERENCIAL (limma) — OPCIONAL, MAS MUITO ÚTIL
# ---------------------------------------------------------------------------
def encontrar_rscript() -> str | None:
    p = shutil.which("Rscript")
    if p:
        return p
    candidatos = sorted(Path("C:/Program Files/R").glob("R-*/bin/Rscript.exe"), reverse=True)
    return str(candidatos[0]) if candidatos else None


def carregar_deg() -> pd.DataFrame:
    csv = TMP / "deg_full.csv"
    if not csv.exists():
        rds = ROOT / "data" / "processed" / "resultado_limma_DCM_vs_CTRL.rds"
        if not rds.exists():
            log("Resultado do limma não encontrado; seguindo sem log2FC.")
            return pd.DataFrame(columns=["gene", "log2FoldChange", "padj", "significant", "direction"])
        rscript = encontrar_rscript()
        if not rscript:
            log("Rscript não encontrado; seguindo sem log2FC.")
            return pd.DataFrame(columns=["gene", "log2FoldChange", "padj", "significant", "direction"])
        cmd = [
            rscript, "-e",
            f'x <- readRDS("{rds.as_posix()}"); write.csv(x$deg, "{csv.as_posix()}", row.names = FALSE)',
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    deg = pd.read_csv(csv)
    return deg[["gene", "log2FoldChange", "padj", "significant", "direction"]]


# ---------------------------------------------------------------------------
# 4) KEGG — VIAS E GENES
# ---------------------------------------------------------------------------
def _strip_kegg_species(nome: str) -> str:
    return re.sub(r"\s*-\s*Homo sapiens \(human\)\s*$", "", nome).strip()


def _parse_kegg_entry(texto: str) -> dict:
    info = {"name": None, "class": None, "entry": None}
    for linha in texto.splitlines():
        if linha.startswith("NAME"):
            info["name"] = _strip_kegg_species(linha[12:].strip())
        elif linha.startswith("CLASS"):
            info["class"] = linha[12:].strip()
        elif linha.startswith("ENTRY"):
            info["entry"] = linha[12:].strip()
    return info


def baixar_kegg() -> tuple[pd.DataFrame, dict[str, set[str]], dict[str, str]]:
    log("Consultando KEGG...")
    # Vias que contêm NQO1.
    link_nqo1 = cached_text(
        "kegg_nqo1_link_pathway.tsv",
        "https://rest.kegg.jp/link/pathway/hsa:1728",
    )
    pathway_ids = []
    for linha in link_nqo1.strip().splitlines():
        partes = linha.split("\t")
        if len(partes) == 2:
            pathway_ids.append(partes[1].replace("path:", "").strip())
    pathway_ids = list(dict.fromkeys(pathway_ids))

    # Lista completa de genes humanos do KEGG: ID -> símbolo e símbolo -> ID.
    listagem = cached_text("kegg_list_hsa.tsv", "https://rest.kegg.jp/list/hsa")
    id2sym: dict[str, str] = {}
    sym2id: dict[str, str] = {}
    for linha in listagem.strip().splitlines():
        campos = linha.split("\t")
        if len(campos) < 2:
            continue
        kid, desc = campos[0], campos[-1]
        if ";" not in desc:
            continue
        simbolos = [s.strip() for s in desc.split(";")[0].split(",") if s.strip()]
        if not simbolos:
            continue
        id2sym[kid] = simbolos[0]
        for s in simbolos:
            sym2id.setdefault(s, kid)

    # Genes de cada via.
    genes_por_via: dict[str, set[str]] = {}
    for pid in pathway_ids:
        texto = cached_text(
            f"kegg_link_hsa_{pid}.tsv",
            f"https://rest.kegg.jp/link/hsa/{pid}",
        )
        genes = set()
        for linha in texto.strip().splitlines():
            partes = linha.split("\t")
            if len(partes) != 2:
                continue
            kid = partes[1].strip()
            if kid in id2sym:
                genes.add(id2sym[kid])
        genes_por_via[pid] = genes

    # Nome e classe de cada via.
    linhas = []
    for pid in pathway_ids:
        texto = cached_text(f"kegg_get_{pid}.txt", f"https://rest.kegg.jp/get/{pid}")
        info = _parse_kegg_entry(texto)
        classe = info["class"]
        if not classe:
            # Vias globais/overview não trazem a linha CLASS no KEGG REST.
            if pid == "hsa01100":
                classe = "Metabolism; Global and overview maps"
            elif pid == "hsa01240":
                classe = "Metabolism; Metabolism of cofactors and vitamins"
            else:
                classe = "Não classificada"
        linhas.append(
            {
                "pathway_id": pid,
                "pathway_name": info["name"] or pid,
                "class": classe,
                "class_group": classe.split(";")[0].strip(),
                "n_genes_pathway": len(genes_por_via.get(pid, set())),
                "url": f"https://www.kegg.jp/pathway/{pid}",
            }
        )
    kegg = pd.DataFrame(linhas)
    log(f"KEGG: {len(kegg)} vias; classes: {', '.join(sorted(kegg['class_group'].unique()))}.")
    return kegg, genes_por_via, sym2id


# ---------------------------------------------------------------------------
# 5) REACTOME — VIAS, HIERARQUIA E PARTICIPANTES
# ---------------------------------------------------------------------------
def _extrair_genes_reactome(participantes: list) -> set[str]:
    genes = set()
    for p in participantes:
        for ref in p.get("refEntities", []) or []:
            dn = ref.get("displayName", "") or ""
            # Ex.: "UniProt:P15559 NQO1" ou "ENSEMBL:ENSG00000181019 NQO1"
            m = re.match(r"^(?:UniProt|ENSEMBL):\S+\s+(.+)$", dn)
            if m:
                genes.add(m.group(1).strip())
    return genes


def _cadeia_ancestrais(ancestrais: list) -> list[str]:
    if not ancestrais:
        return []
    # A API devolve uma lista de caminhos; usamos o primeiro.
    caminho = ancestrais[0] if isinstance(ancestrais[0], list) else ancestrais
    return [x.get("displayName", "") for x in caminho if x.get("displayName")]


def baixar_reactome() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, set[str]], str]:
    log("Consultando Reactome...")
    busca = cached_json(
        "reactome_search_nqo1.json",
        "https://reactome.org/ContentService/search/query",
        params={"query": GENE, "species": "Homo sapiens", "cluster": "true"},
    )
    db_id = None
    for resultado in busca.get("results", []):
        for entrada in resultado.get("entries", []):
            if (
                entrada.get("type") == "Protein"
                and entrada.get("referenceIdentifier") == UNIPROT
            ):
                db_id = str(entrada["dbId"])
                break
        if db_id:
            break
    if not db_id:
        raise RuntimeError("Não consegui mapear NQO1 para uma entidade do Reactome.")

    vias_baixas = cached_json(
        "reactome_nqo1_low_pathways.json",
        f"https://reactome.org/ContentService/data/pathways/low/entity/{db_id}/allForms",
        params={"species": SPECIES},
    )

    genes_por_via: dict[str, set[str]] = {}
    linhas = []
    hierarquia = []
    for via in vias_baixas:
        vid = str(via["dbId"])
        stid = via["stId"]
        nome = via.get("displayName") or via.get("name", [stid])[0]

        ancestrais = cached_json(
            f"reactome_ancestors_{vid}.json",
            f"https://reactome.org/ContentService/data/event/{vid}/ancestors",
        )
        cadeia = _cadeia_ancestrais(ancestrais)
        if not cadeia:
            cadeia = [nome]
        top_level = cadeia[-1]

        participantes = cached_json(
            f"reactome_participants_{vid}.json",
            f"https://reactome.org/ContentService/data/participants/{vid}",
        )
        genes = _extrair_genes_reactome(participantes)
        genes_por_via[stid] = genes

        linhas.append(
            {
                "pathway_id": stid,
                "pathway_name": nome,
                "top_level": top_level,
                "ancestor_chain": " > ".join(cadeia),
                "n_participants": len(genes),
                "n_ancestors": len(cadeia),
                "url": f"https://reactome.org/content/detail/{stid}",
            }
        )

        for ordem, evento in enumerate(cadeia):
            hierarquia.append(
                {
                    "pathway_id": stid,
                    "pathway_name": nome,
                    "ordem": ordem,
                    "evento": evento,
                    "direto": ordem == 0,
                    "top_level": top_level,
                    "url": f"https://reactome.org/content/detail/{stid}",
                }
            )

    reactome = pd.DataFrame(linhas)
    hierarquia_df = pd.DataFrame(hierarquia)
    log(f"Reactome: {len(reactome)} vias de baixo nível; "
        f"top-level: {', '.join(sorted(reactome['top_level'].unique()))}.")
    return reactome, hierarquia_df, genes_por_via, db_id


# ---------------------------------------------------------------------------
# 6) CRUZAMENTOS E TABELAS
# ---------------------------------------------------------------------------
def anotar_nos(
    g: nx.Graph,
    arestas: pd.DataFrame,
    kegg: pd.DataFrame,
    genes_kegg: dict[str, set[str]],
    reactome: pd.DataFrame,
    genes_reactome: dict[str, set[str]],
    deg: pd.DataFrame,
) -> pd.DataFrame:
    # Score direto NQO1-parceiro.
    diretos = arestas[arestas["is_nqo1_edge"]].copy()
    score_direto = {}
    for _, r in diretos.iterrows():
        parceiro = r["gene_b"] if r["gene_a"] == GENE else r["gene_a"]
        score_direto[parceiro] = float(r["score"])
    # O próprio NQO1 recebe o maior score observado entre suas arestas diretas.
    score_direto[GENE] = float(diretos["score"].max()) if not diretos.empty else np.nan

    entre = nx.betweenness_centrality(g, weight=None)
    prox = nx.closeness_centrality(g)
    agrup = nx.clustering(g)

    # ID STRING de cada nó (o nó pode aparecer como gene_a ou gene_b).
    string_ids: dict[str, str] = {}
    for _, r in arestas.iterrows():
        string_ids[r["gene_a"]] = r["string_id_a"]
        string_ids[r["gene_b"]] = r["string_id_b"]

    linhas = []
    for no in g.nodes:
        match = ALIAS.get(no, no)
        vias_kegg = [
            pid for pid, genes in genes_kegg.items()
            if match in genes
        ]
        vias_reactome = [
            pid for pid, genes in genes_reactome.items()
            if match in genes
        ]
        linhas.append(
            {
                "gene": no,
                "is_nqo1": no == GENE,
                "string_id": string_ids.get(no),
                "degree": g.degree(no),
                "weighted_degree": round(sum(d["score"] for _, _, d in g.edges(no, data=True)), 3),
                "direct_score_nqo1": score_direto.get(no, np.nan),
                "betweenness": round(entre.get(no, np.nan), 4),
                "closeness": round(prox.get(no, np.nan), 4),
                "clustering": round(agrup.get(no, np.nan), 4),
                "n_kegg": len(vias_kegg),
                "n_reactome": len(vias_reactome),
                "n_pathways_total": len(vias_kegg) + len(vias_reactome),
                "kegg_pathway_ids": ";".join(sorted(vias_kegg)),
                "reactome_pathway_ids": ";".join(sorted(vias_reactome)),
            }
        )
    nodes = pd.DataFrame(linhas)

    if not deg.empty:
        nodes = nodes.merge(deg, on="gene", how="left")
    else:
        for c in ["log2FoldChange", "padj", "significant", "direction"]:
            nodes[c] = np.nan

    # Nome legível das vias (para o relatório), sem repetir a tabela inteira.
    kegg_nome = dict(zip(kegg["pathway_id"], kegg["pathway_name"]))
    react_nome = dict(zip(reactome["pathway_id"], reactome["pathway_name"]))
    nodes["kegg_pathway_names"] = nodes["kegg_pathway_ids"].apply(
        lambda s: "; ".join(kegg_nome.get(x, x) for x in str(s).split(";") if x)
    )
    nodes["reactome_pathway_names"] = nodes["reactome_pathway_ids"].apply(
        lambda s: "; ".join(react_nome.get(x, x) for x in str(s).split(";") if x)
    )

    # Genes do dataset sem medida ficam NaN; para desenhar usamos cinza.
    nodes["medido"] = nodes["log2FoldChange"].notna()
    nodes = nodes.sort_values(["is_nqo1", "degree"], ascending=[False, False]).reset_index(drop=True)
    return nodes


def completar_kegg_reactome(
    kegg: pd.DataFrame,
    genes_kegg: dict[str, set[str]],
    reactome: pd.DataFrame,
    genes_reactome: dict[str, set[str]],
    parceiros: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    def lista_parceiros(genes: set[str]) -> list[str]:
        return sorted(
            p for p in parceiros
            if ALIAS.get(p, p) in genes
        )

    kegg = kegg.copy()
    kegg["string_partners"] = kegg["pathway_id"].map(
        lambda pid: ";".join(lista_parceiros(genes_kegg.get(pid, set())))
    )
    kegg["n_string_partners"] = kegg["string_partners"].apply(
        lambda s: 0 if not s else len(s.split(";"))
    )
    kegg = kegg.sort_values(
        ["n_string_partners", "pathway_name"], ascending=[False, True]
    ).reset_index(drop=True)

    reactome = reactome.copy()
    reactome["string_partners"] = reactome["pathway_id"].map(
        lambda pid: ";".join(lista_parceiros(genes_reactome.get(pid, set())))
    )
    reactome["n_string_partners"] = reactome["string_partners"].apply(
        lambda s: 0 if not s else len(s.split(";"))
    )
    reactome = reactome.sort_values(
        ["n_string_partners", "pathway_name"], ascending=[False, True]
    ).reset_index(drop=True)
    return kegg, reactome


# ---------------------------------------------------------------------------
# 7) FIGURAS
# ---------------------------------------------------------------------------
def _posicoes_radiais(g: nx.Graph, nodes: pd.DataFrame) -> dict:
    parceiros = [n for n in g.nodes if n != GENE]
    ordem = nodes[nodes["gene"] != GENE].sort_values(
        ["degree", "gene"], ascending=[False, True]
    )["gene"].tolist()
    # Garante que todos os nós do grafo estejam na ordem.
    parceiros = [n for n in ordem if n in parceiros] + [n for n in parceiros if n not in ordem]
    angulos = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, len(parceiros), endpoint=False)
    pos = {GENE: np.array([0.0, 0.0])}
    for ang, no in zip(angulos, parceiros):
        pos[no] = np.array([np.cos(ang), np.sin(ang)])
    return pos


def _fonte(tamanho: int):
    for nome in ("arial.ttf", "DejaVuSans.ttf", "LiberationSans-Regular.ttf"):
        try:
            return ImageFont.truetype(nome, tamanho)
        except Exception:
            continue
    return ImageFont.load_default()


def figura_rede(g: nx.Graph, nodes: pd.DataFrame, arestas: pd.DataFrame) -> None:
    log("Desenhando a rede PPI...")
    pos = _posicoes_radiais(g, nodes)
    info = nodes.set_index("gene")

    fig, ax = plt.subplots(figsize=(16, 14))

    # Arestas entre parceiros (fundinho, para mostrar módulos sem poluir).
    seg_parceiros = []
    for a, b, d in g.edges(data=True):
        if GENE not in (a, b):
            seg_parceiros.append((pos[a], pos[b]))
    if seg_parceiros:
        ax.add_collection(
            LineCollection(
                seg_parceiros,
                colors=[(0.55, 0.55, 0.55, 0.075)],
                linewidths=0.45,
                zorder=1,
            )
        )

    # Arestas diretas NQO1-parceiro.
    seg_diretas = []
    larguras = []
    for a, b, d in g.edges(data=True):
        if GENE in (a, b):
            seg_diretas.append((pos[a], pos[b]))
            # Score 0.70 -> 0.6 pt; score 1.00 -> 2.8 pt.
            larguras.append(0.6 + 2.2 * max(0.0, (d["score"] - 0.70) / 0.30))
    ax.add_collection(
        LineCollection(
            seg_diretas,
            colors=[(0.28, 0.28, 0.28, 0.55)],
            linewidths=larguras,
            zorder=2,
        )
    )

    # Nós.
    valores = info.loc[list(g.nodes), "log2FoldChange"].astype(float).values
    medido = ~np.isnan(valores)
    lim_fc = np.nanmax(np.abs(valores)) if medido.any() else 1.0
    lim_fc = max(1.0, np.ceil(lim_fc * 10) / 10)
    norm = Normalize(vmin=-lim_fc, vmax=lim_fc)
    cmap = plt.get_cmap("RdBu_r")

    graus = np.array([g.degree(n) for n in g.nodes])
    tamanhos = np.array([2300 if n == GENE else 120 + 42 * g.degree(n) for n in g.nodes])

    # Nós sem medida em cinza; nós medidos com gradiente.
    idx_cinza = [i for i, n in enumerate(g.nodes) if np.isnan(info.loc[n, "log2FoldChange"])]
    idx_cor = [i for i, n in enumerate(g.nodes) if not np.isnan(info.loc[n, "log2FoldChange"])]
    xs = [pos[n][0] for n in g.nodes]
    ys = [pos[n][1] for n in g.nodes]

    if idx_cinza:
        ax.scatter(
            [xs[i] for i in idx_cinza], [ys[i] for i in idx_cinza],
            s=[tamanhos[i] for i in idx_cinza],
            c="#c9c9c9", edgecolors="#777777", linewidths=0.6, zorder=3,
        )
    if idx_cor:
        sc = ax.scatter(
            [xs[i] for i in idx_cor], [ys[i] for i in idx_cor],
            s=[tamanhos[i] for i in idx_cor],
            c=[valores[i] for i in idx_cor],
            cmap=cmap, norm=norm,
            edgecolors="#222222", linewidths=0.7, zorder=4,
        )
    else:
        sc = None

    # Contorno especial do NQO1.
    ax.scatter(
        [0], [0], s=[2600], facecolors="none",
        edgecolors="black", linewidths=2.6, zorder=5,
    )

    # Rótulos externos.
    top_hubs = set(nodes[nodes["gene"] != GENE].nlargest(8, "degree")["gene"])
    for no in g.nodes:
        if no == GENE:
            continue
        x, y = pos[no]
        ha = "left" if x > 0.05 else ("right" if x < -0.05 else "center")
        peso = "bold" if no in top_hubs else "normal"
        ax.text(
            x * 1.105, y * 1.105, no,
            ha=ha, va="center", fontsize=8.4, fontweight=peso,
            color="#111111",
            path_effects=[pe.withStroke(linewidth=2.6, foreground="white")],
            zorder=6,
        )
    ax.text(
        0, 0, GENE, ha="center", va="center",
        fontsize=11, fontweight="bold", color="white", zorder=7,
    )

    # Legenda de tamanho.
    handles_tam = [
        Line2D([0], [0], marker="o", color="none", label=f"grau {d}",
               markerfacecolor="#cccccc", markeredgecolor="#555555",
               markersize=np.sqrt(120 + 42 * d) / 2.2)
        for d in (5, 15, 25)
    ]
    leg1 = ax.legend(
        handles=handles_tam, loc="lower left", frameon=True, framealpha=0.9,
        title="Tamanho do nó", fontsize=8, title_fontsize=9, borderpad=0.7,
    )
    ax.add_artist(leg1)

    # Legenda de arestas.
    handles_arestas = [
        Line2D([0], [0], color=(0.28, 0.28, 0.28, 0.7), lw=2.2,
               label="NQO1–parceiro (score alto)"),
        Line2D([0], [0], color=(0.55, 0.55, 0.55, 0.35), lw=1.2,
               label="parceiro–parceiro"),
    ]
    ax.legend(
        handles=handles_arestas, loc="lower right", frameon=True, framealpha=0.9,
        title="Arestas STRING", fontsize=8, title_fontsize=9, borderpad=0.7,
    )

    if sc is not None:
        cax = fig.add_axes([0.435, 0.055, 0.18, 0.018])
        cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
        cb.set_label("log2FC (DCM vs CTRL)", fontsize=9)
        cb.ax.tick_params(labelsize=8)

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.32, 1.32)
    ax.set_ylim(-1.25, 1.25)
    ax.set_title(
        "Rede PPI centrada em NQO1 no STRING",
        fontsize=18, fontweight="bold", pad=18,
    )
    ax.text(
        0.5, 0.965,
        f"Homo sapiens | score ≥ {STRING_SCORE} | {g.number_of_nodes()} nós e "
        f"{g.number_of_edges()} arestas | {g.degree(GENE)} parceiros diretos",
        transform=ax.transAxes, ha="center", va="top", fontsize=10.5, color="#444444",
    )
    fig.text(
        0.5, 0.018,
        "Anel externo: parceiros diretos de NQO1. Cor: log2FC na DCM (azul = down, vermelho = up; "
        "cinza = gene não medido no dataset).\n"
        "Arestas escuras: interação NQO1–parceiro; arestas claras: interação entre parceiros. "
        "Tamanho do nó: grau na rede.",
        ha="center", va="bottom", fontsize=8.2, color="#555555",
    )
    fig.subplots_adjust(top=0.93, bottom=0.075, left=0.02, right=0.98)
    fig.savefig(FIG / "nqo1_rede_ppi_string.png", dpi=300)
    plt.close(fig)
    log("Figura da rede salva.")


def _quebrar(texto: str, largura: int = 42) -> str:
    return "\n".join(
        texto[i:i + largura] for i in range(0, len(texto), largura)
    )


def figura_kegg(kegg: pd.DataFrame) -> None:
    log("Desenhando barras KEGG...")
    df = kegg.sort_values(["n_string_partners", "pathway_name"], ascending=[True, True]).copy()
    df["rotulo"] = [
        f"{_quebrar(nome, 48)}\n[{pid}]"
        for nome, pid in zip(df["pathway_name"], df["pathway_id"])
    ]
    cores = {
        "Metabolism": "#2A9D8F",
        "Human Diseases": "#E76F51",
        "Não classificada": "#888888",
    }
    df["cor"] = df["class_group"].map(cores).fillna("#888888")

    fig, ax = plt.subplots(figsize=(13, 7.2))
    barras = ax.barh(
        df["rotulo"], df["n_string_partners"],
        color=df["cor"], edgecolor="#333333", linewidth=0.5, height=0.68,
    )
    for barra, (_, r) in zip(barras, df.iterrows()):
        x = barra.get_width()
        ax.text(
            x + 0.15, barra.get_y() + barra.get_height() / 2,
            f"{int(r['n_string_partners'])} parceiros  "
            f"({int(r['n_genes_pathway'])} genes na via)",
            va="center", ha="left", fontsize=8.5, color="#333333",
        )
    ax.set_xlabel("Nº de parceiros diretos de NQO1 (STRING) presentes na via", fontsize=10)
    ax.set_xlim(0, max(df["n_string_partners"].max() + 3.2, 3))
    ax.set_title(
        "Vias KEGG que contêm NQO1",
        fontsize=17, fontweight="bold", pad=12,
    )
    ax.text(
        0.0, 1.01,
        "Barras = quantos parceiros da rede STRING também estão anotados na via. "
        "Entre parênteses: total de genes da via no KEGG.",
        transform=ax.transAxes, fontsize=9, color="#555555", va="bottom",
    )
    leg = [
        Patch(facecolor=cores["Metabolism"], edgecolor="#333333", label="Metabolism"),
        Patch(facecolor=cores["Human Diseases"], edgecolor="#333333", label="Human Diseases"),
    ]
    ax.legend(handles=leg, loc="lower right", frameon=True, fontsize=8.5, title="Classe KEGG")
    ax.tick_params(axis="y", labelsize=8.3)
    ax.tick_params(axis="x", labelsize=9)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "nqo1_kegg_pathways.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log("Figura KEGG salva.")


def figura_reactome(reactome: pd.DataFrame) -> None:
    log("Desenhando barras Reactome...")
    df = reactome.sort_values(["n_string_partners", "pathway_name"], ascending=[True, True]).copy()
    df["rotulo"] = [
        f"{_quebrar(nome, 44)}\n[{pid}]"
        for nome, pid in zip(df["pathway_name"], df["pathway_id"])
    ]
    palette = {
        "Metabolism": "#264653",
        "Chromatin organization": "#8E7DBE",
        "Cellular responses to stimuli": "#E9C46A",
        "Cellular responses to stress": "#F4A261",
        "Cellular response to chemical stress": "#E76F51",
    }
    df["cor"] = df["top_level"].map(palette).fillna("#888888")

    fig, ax = plt.subplots(figsize=(13, 5.4))
    barras = ax.barh(
        df["rotulo"], df["n_string_partners"],
        color=df["cor"], edgecolor="#333333", linewidth=0.5, height=0.6,
    )
    for barra, (_, r) in zip(barras, df.iterrows()):
        x = barra.get_width()
        ax.text(
            x + 0.12, barra.get_y() + barra.get_height() / 2,
            f"{int(r['n_string_partners'])} parceiros  ({int(r['n_participants'])} genes/participantes)",
            va="center", ha="left", fontsize=8.8, color="#333333",
        )
    ax.set_xlabel("Nº de parceiros diretos de NQO1 (STRING) também anotados na via", fontsize=10)
    ax.set_xlim(0, max(df["n_string_partners"].max() + 3.0, 3))
    ax.set_title(
        "Vias Reactome que contêm NQO1",
        fontsize=17, fontweight="bold", pad=12,
    )
    ax.text(
        0.0, 1.015,
        "Barras = parceiros STRING também presentes na via. "
        "Cor = via de nível superior (top-level) do Reactome.",
        transform=ax.transAxes, fontsize=9, color="#555555", va="bottom",
    )
    handles = [
        Patch(facecolor=palette.get(x, "#888888"), edgecolor="#333333", label=x)
        for x in df["top_level"].unique()
    ]
    ax.legend(handles=handles, loc="lower right", frameon=True, fontsize=8.2,
              title="Top-level Reactome")
    ax.tick_params(axis="y", labelsize=8.5)
    ax.tick_params(axis="x", labelsize=9)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "nqo1_reactome_pathways.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log("Figura Reactome salva.")


def figura_heatmap(nodes: pd.DataFrame, kegg: pd.DataFrame, reactome: pd.DataFrame) -> None:
    log("Desenhando heatmap de vias...")
    # Só entram nós que compartilham pelo menos uma via com NQO1.
    df = nodes[(nodes["is_nqo1"]) | (nodes["n_pathways_total"] > 0)].copy()
    df = df.sort_values(["is_nqo1", "n_pathways_total", "degree"], ascending=[False, False, False])
    genes = df["gene"].tolist()

    colunas = []
    grupos = []
    for _, r in kegg.sort_values("pathway_id").iterrows():
        nome_curto = _quebrar(r["pathway_name"], 26).replace("\n", " ")
        if len(nome_curto) > 26:
            nome_curto = nome_curto[:24] + "…"
        colunas.append(f"{r['pathway_id']}\n{nome_curto}")
        grupos.append("KEGG")
    for _, r in reactome.sort_values("pathway_id").iterrows():
        nome_curto = _quebrar(r["pathway_name"], 26).replace("\n", " ")
        if len(nome_curto) > 26:
            nome_curto = nome_curto[:24] + "…"
        colunas.append(f"{r['pathway_id']}\n{nome_curto}")
        grupos.append("Reactome")

    ids_kegg = kegg.sort_values("pathway_id")["pathway_id"].tolist()
    ids_react = reactome.sort_values("pathway_id")["pathway_id"].tolist()

    matriz = np.zeros((len(genes), len(ids_kegg) + len(ids_react)), dtype=float)
    for i, gene in enumerate(genes):
        k_ids = set(str(nodes.loc[nodes["gene"] == gene, "kegg_pathway_ids"].iloc[0]).split(";"))
        r_ids = set(str(nodes.loc[nodes["gene"] == gene, "reactome_pathway_ids"].iloc[0]).split(";"))
        for j, pid in enumerate(ids_kegg):
            matriz[i, j] = 1.0 if pid in k_ids else 0.0
        for j, pid in enumerate(ids_react):
            matriz[i, len(ids_kegg) + j] = 1.0 if pid in r_ids else 0.0

    altura = max(5.5, 0.28 * len(genes) + 2.2)
    fig, ax = plt.subplots(figsize=(13.5, altura))
    cmap = matplotlib.colors.ListedColormap(["#f7f7f7", "#1f6f8b"])
    ax.imshow(matriz, aspect="auto", cmap=cmap, vmin=0, vmax=1)

    ax.set_xticks(range(len(colunas)))
    ax.set_xticklabels(colunas, rotation=55, ha="right", fontsize=7.2)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8.2)
    for tick, gene in zip(ax.get_yticklabels(), genes):
        if gene == GENE:
            tick.set_fontweight("bold")

    # Marcações e grade.
    for i in range(len(genes)):
        for j in range(len(colunas)):
            if matriz[i, j] > 0:
                ax.text(j, i, "•", ha="center", va="center", color="white", fontsize=8)
    ax.set_xticks(np.arange(-0.5, len(colunas), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(genes), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Separador KEGG / Reactome.
    if ids_kegg and ids_react:
        ax.axvline(len(ids_kegg) - 0.5, color="#333333", linewidth=2)
        ax.text(len(ids_kegg) / 2 - 0.5, -1.1, "KEGG", ha="center", va="bottom",
                fontsize=10, fontweight="bold", color="#2A9D8F")
        ax.text(len(ids_kegg) + len(ids_react) / 2 - 0.5, -1.1, "Reactome",
                ha="center", va="bottom", fontsize=10, fontweight="bold", color="#8E7DBE")

    ax.set_title(
        "Parceiros de NQO1 que compartilham vias KEGG/Reactome",
        fontsize=15, fontweight="bold", pad=30,
    )
    ax.set_xlabel("Via (ID e nome abreviado)", fontsize=9, labelpad=8)
    ax.set_ylabel("Proteína", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "nqo1_heatmap_vias_parceiros.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log("Heatmap salvo.")


def figura_kegg_diagrama() -> None:
    log("Baixando e destacando o mapa KEGG hsa05208...")
    try:
        imagem_bytes = cached_bytes(
            "kegg_hsa05208.png",
            "https://rest.kegg.jp/get/hsa05208/image",
        )
        html = cached_text(
            "kegg_hsa05208_highlight.html",
            "https://www.kegg.jp/kegg-bin/show_pathway?hsa05208+1728",
        )
        # Encontra as áreas do NQO1 no mapa.
        coords = []
        for tag in re.findall(r"<area\b[^>]*>", html):
            if "NQO1" in tag and "coords=" in tag:
                m = re.search(r'coords="([^"]+)"', tag)
                if m:
                    partes = [int(x) for x in re.findall(r"\d+", m.group(1))]
                    if len(partes) >= 4:
                        coords.append(tuple(partes[:4]))
        if not coords:
            log("Não encontrei as coordenadas do NQO1 no HTML do KEGG; salvando mapa sem destaque.")
            (FIG / "nqo1_kegg_diagram_hsa05208.png").write_bytes(imagem_bytes)
            return

        img = Image.open(io.BytesIO(imagem_bytes)).convert("RGB")
        draw = ImageDraw.Draw(img)
        fonte = _fonte(30)
        for i, (x1, y1, x2, y2) in enumerate(coords):
            draw.rectangle([x1 - 7, y1 - 7, x2 + 7, y2 + 7], outline=(220, 20, 20), width=5)
            if i == 0:
                draw.text((x1 - 7, max(0, y1 - 42)), "NQO1", fill=(220, 20, 20), font=fonte)
        img.save(FIG / "nqo1_kegg_diagram_hsa05208.png")
        log(f"Mapa KEGG salvo com {len(coords)} destaque(s) para NQO1.")
    except Exception as e:
        log(f"Não foi possível gerar o mapa KEGG: {e}")


def figura_reactome_diagrama(reactome: pd.DataFrame) -> None:
    log("Baixando o diagrama Reactome da via NFE2L2...")
    try:
        # A via R-HSA-9818027 tem diagrama; as outras podem não ter.
        candidata = reactome[reactome["pathway_id"] == "R-HSA-9818027"]
        if candidata.empty:
            return
        imagem = cached_bytes(
            "reactome_R-HSA-9818027.png",
            "https://reactome.org/ContentService/exporter/diagram/R-HSA-9818027.png",
        )
        (FIG / "nqo1_reactome_diagram_nfe2l2.png").write_bytes(imagem)
        log("Diagrama Reactome salvo.")
    except Exception as e:
        log(f"Não foi possível baixar o diagrama Reactome: {e}")


# ---------------------------------------------------------------------------
# 8) RELATÓRIO EM MARKDOWN
# ---------------------------------------------------------------------------
def _md_tabela(df: pd.DataFrame, colunas: list[str], renomear: dict[str, str] | None = None) -> str:
    d = df[colunas].copy()
    if renomear:
        d = d.rename(columns=renomear)
    for c in d.columns:
        if d[c].dtype == object:
            d[c] = d[c].fillna("").astype(str)
    return d.to_markdown(index=False)


def escrever_relatorio(
    nodes: pd.DataFrame,
    kegg: pd.DataFrame,
    reactome: pd.DataFrame,
    hierarquia: pd.DataFrame,
    arestas: pd.DataFrame,
    g: nx.Graph,
) -> None:
    log("Escrevendo relatório Markdown...")
    nqo = nodes[nodes["is_nqo1"]].iloc[0]
    top = nodes[nodes["gene"] != GENE].nlargest(15, "degree")

    top_fmt = top.copy()
    top_fmt["direct_score_nqo1"] = top_fmt["direct_score_nqo1"].map(
        lambda x: "—" if pd.isna(x) else f"{x:.3f}"
    )
    top_fmt["log2FoldChange"] = top_fmt["log2FoldChange"].map(
        lambda x: "—" if pd.isna(x) else f"{x:+.2f}"
    )
    top_fmt["padj"] = top_fmt["padj"].map(
        lambda x: "—" if pd.isna(x) else f"{x:.2e}"
    )
    top_md = _md_tabela(
        top_fmt,
        ["gene", "degree", "direct_score_nqo1", "log2FoldChange", "padj", "n_kegg", "n_reactome"],
        {
            "gene": "Gene",
            "degree": "Grau",
            "direct_score_nqo1": "Score NQO1",
            "log2FoldChange": "log2FC",
            "padj": "padj",
            "n_kegg": "Vias KEGG*",
            "n_reactome": "Vias Reactome*",
        },
    )

    kegg_md = _md_tabela(
        kegg,
        ["pathway_id", "pathway_name", "class", "n_genes_pathway", "n_string_partners", "string_partners"],
        {
            "pathway_id": "ID",
            "pathway_name": "Via KEGG",
            "class": "Classe",
            "n_genes_pathway": "Genes na via",
            "n_string_partners": "Parceiros STRING",
            "string_partners": "Parceiros",
        },
    )
    react_md = _md_tabela(
        reactome,
        ["pathway_id", "pathway_name", "top_level", "n_participants", "n_string_partners", "string_partners"],
        {
            "pathway_id": "ID",
            "pathway_name": "Via Reactome",
            "top_level": "Top-level",
            "n_participants": "Genes/participantes",
            "n_string_partners": "Parceiros STRING",
            "string_partners": "Parceiros",
        },
    )

    hier_md = _md_tabela(
        hierarquia[hierarquia["direto"] == False],
        ["pathway_id", "evento", "top_level"],
        {
            "pathway_id": "Via NQO1",
            "evento": "Evento ancestral",
            "top_level": "Top-level",
        },
    )

    metricas = pd.DataFrame(
        [
            ("Gene alvo", GENE),
            ("UniProt", UNIPROT),
            ("Espécie", "Homo sapiens (9606)"),
            ("Rede STRING — nós", g.number_of_nodes()),
            ("Rede STRING — arestas", g.number_of_edges()),
            ("Parceiros diretos de NQO1", g.degree(GENE)),
            ("Score mínimo STRING", STRING_SCORE),
            ("Vias KEGG com NQO1", len(kegg)),
            ("Vias Reactome (baixo nível) com NQO1", len(reactome)),
            ("NQO1 log2FC (DCM vs CTRL)", f"{nqo['log2FoldChange']:.3f}"),
            ("NQO1 padj", f"{nqo['padj']:.2e}"),
            ("NQO1 direção", nqo["direction"]),
        ],
        columns=["Métrica", "Valor"],
    )
    metricas.to_csv(TAB / "nqo1_resumo_metricas.csv", index=False)
    metricas_md = metricas.to_markdown(index=False)

    texto = f"""# NQO1 — rede PPI no STRING e vias no KEGG e Reactome

**Gene alvo:** NQO1 (UniProt {UNIPROT}) · **Espécie:** *Homo sapiens* (9606)  
**Dataset de expressão:** GSE116250 (DCM vs CTRL, limma)  
**Data da análise:** {pd.Timestamp.today().strftime('%d/%m/%Y')}

---

## 1. Resumo

| Banco | O que foi consultado | Resultado principal |
|---|---|---|
| **STRING** | Rede de vizinhança direta de NQO1, score ≥ {STRING_SCORE} | **{g.number_of_nodes()} nós** e **{g.number_of_edges()} arestas**; NQO1 tem **{g.degree(GENE)} parceiros diretos** |
| **KEGG** | Vias que contêm NQO1 e genes dessas vias | **{len(kegg)} vias**; NQO1 aparece em vias de metabolismo e de doenças humanas |
| **Reactome** | Vias de baixo nível que contêm NQO1 e hierarquia | **{len(reactome)} vias**; a principal é a via **NFE2L2/KEAP1** |

> **NQO1 no dataset:** log2FC = **{nqo['log2FoldChange']:.2f}** (padj = **{nqo['padj']:.2e}**), ou seja, **{nqo['direction']}regulado** na cardiomiopatia dilatada (DCM) em relação ao controle.

![Rede PPI centrada em NQO1](output/figures/nqo1_rede_ppi_string.png)

**Como ler a rede:** cada bolinha é uma proteína; o centro é o NQO1. O anel externo mostra os {g.degree(GENE)} parceiros diretos. O tamanho do nó representa o grau (nº de conexões na rede); a cor representa o log2FC na DCM (azul = down, vermelho = up; cinza = gene não medido no dataset). As linhas escuras ligam NQO1 aos parceiros; as linhas claras mostram interações entre os próprios parceiros.

---

## 2. Métodos (resumo)

1. **STRING** — API `network` (`species=9606`, `required_score={STRING_SCORE}`, `limit={STRING_LIMIT}`). O score combinado varia de 0 a 1; ≥ 0,700 é considerado alta confiança. A rede inclui NQO1, seus parceiros diretos e as interações entre esses parceiros.
2. **KEGG** — API REST (`/link/pathway/hsa:1728`, `/link/hsa/{{pathway}}`, `/get/{{pathway}}`). Foram obtidas as vias que contêm NQO1, seus nomes, classes e genes.
3. **Reactome** — Content Service (`/search/query`, `/data/pathways/low/entity/...`, `/data/event/.../ancestors`, `/data/participants/...`). Foram obtidas as vias de baixo nível, a hierarquia até o top-level e os participantes.
4. **Expressão** — resultado do `limma` do projeto (`resultado_limma_DCM_vs_CTRL.rds`), usado apenas para colorir os nós.
5. **Cruzamento** — cada parceiro do STRING foi comparado com os genes/participantes das vias KEGG e Reactome que contêm NQO1.

> **Importante:** KEGG e Reactome **não são bases de PPI**. A rede de interações proteína-proteína vem do STRING; KEGG/Reactome entram como anotação de vias e contexto funcional.

---

## 3. Resultados

### 3.1 Rede PPI do STRING

A rede tem **{g.number_of_nodes()} proteínas** e **{g.number_of_edges()} interações** de alta confiança. O NQO1 é o nó central, com **{g.degree(GENE)} parceiros diretos**. Os nós mais conectados (hubs) são **{', '.join(top.head(5)['gene'].tolist())}**.

**Tabela 1 — 15 parceiros com maior grau na rede** (o NQO1 não entra nesta lista; ele é o centro).

{top_md}

\\* Número de vias **de NQO1** (KEGG ou Reactome) em que o parceiro também aparece; não é o total de vias do gene.

**Arquivos completos:**
- `output/tables/nqo1_string_ppi_arestas.csv` — todas as arestas, com score e canais de evidência do STRING (pode ser importado no Cytoscape como edge list).
- `output/tables/nqo1_string_ppi_nos.csv` — todos os nós, com grau, centralidades, log2FC e vias compartilhadas.

#### Módulos que aparecem na rede

- **Módulo NRF2/KEAP1:** NFE2L2, KEAP1, GCLC, GCLM, HMOX1, GSR, TXNRD1, SRXN1, TXN. É o eixo clássico de resposta antioxidante e detoxificação.
- **Família GPX:** GPX2, GPX3, GPX5, GPX6, GPX7 e GPX8 formam um grupo muito conectado, ligado ao metabolismo de peróxidos.
- **Eixo TP53/STUB1/PARK7:** TP53, TP73, STUB1 e PARK7 aparecem como conexões de estresse celular e controle de qualidade proteica.
- **Metabolismo de quinonas/vitamina K:** GGCX, VKORC1, VKORC1L1, CRYZ e NDUFS7.

### 3.2 KEGG

O NQO1 está anotado em **{len(kegg)} vias KEGG**:

![Vias KEGG que contêm NQO1](output/figures/nqo1_kegg_pathways.png)

**Tabela 2 — Vias KEGG que contêm NQO1.**

{kegg_md}

**Como ler:** a barra mostra quantos parceiros diretos de NQO1 no STRING também estão naquela via. O número entre parênteses é o total de genes da via no KEGG.

**Destaque:** `hsa05208` (*Chemical carcinogenesis — reactive oxygen species*) é particularmente relevante para o projeto, pois conecta NQO1 ao estresse oxidativo e à sinalização NRF2/KEAP1. O mapa dessa via, com o NQO1 destacado em vermelho, está em:

<img src="output/figures/nqo1_kegg_diagram_hsa05208.png" width="720" alt="Mapa KEGG hsa05208 com NQO1 destacado">

*No mapa, os retângulos vermelhos marcam as ocorrências de NQO1.*

### 3.3 Reactome

O NQO1 está anotado em **{len(reactome)} vias de baixo nível** no Reactome:

![Vias Reactome que contêm NQO1](output/figures/nqo1_reactome_pathways.png)

**Tabela 3 — Vias Reactome de baixo nível que contêm NQO1.**

{react_md}

**Hierarquia até o top-level** (eventos ancestrais, sem repetir a própria via):

{hier_md}

**Destaque:** a via **R-HSA-9818027** (*NFE2L2 regulating anti-oxidant/detoxification enzymes*) é a que mais converge com a rede PPI: **{int(reactome.loc[reactome['pathway_id']=='R-HSA-9818027','n_string_partners'].iloc[0])} parceiros** do STRING também participam dela. Ela sobe na hierarquia até **Cellular responses to stimuli**, passando por *KEAP1-NFE2L2 pathway* e *Cellular response to chemical stress*.

<img src="output/figures/nqo1_reactome_diagram_nfe2l2.png" width="720" alt="Diagrama Reactome da via NFE2L2">

*Diagrama da via R-HSA-9818027 no Reactome.*

### 3.4 Integração: parceiros que compartilham vias com NQO1

![Heatmap de vias compartilhadas](output/figures/nqo1_heatmap_vias_parceiros.png)

O heatmap mostra, para cada parceiro de NQO1, em quais vias KEGG (bloco à esquerda) e Reactome (bloco à direita) ele também aparece. As linhas mais “cheias” são os parceiros que mais convergem funcionalmente com NQO1. O NQO1 aparece em todas as colunas por definição (é o gene consultado).

**Tabela 4 — Parceiros com maior número de vias compartilhadas.**

{_md_tabela(nodes[nodes['gene'] != GENE].nlargest(15, 'n_pathways_total'), ['gene', 'degree', 'n_kegg', 'n_reactome', 'n_pathways_total', 'kegg_pathway_ids', 'reactome_pathway_ids'], {'gene': 'Gene', 'degree': 'Grau', 'n_kegg': 'KEGG', 'n_reactome': 'Reactome', 'n_pathways_total': 'Total', 'kegg_pathway_ids': 'IDs KEGG', 'reactome_pathway_ids': 'IDs Reactome'})}

---

## 4. Métricas do resumo

{metricas_md}

---

## 5. Limitações

- O STRING integra **evidências heterogêneas** (texto, experimentos, bancos de dados, coexpressão, vizinhança genômica, fusão gênica). Uma aresta não significa necessariamente contato físico direto.
- O score ≥ {STRING_SCORE} reduz falsos positivos, mas ainda pode incluir associações funcionais indiretas.
- O KEGG e o Reactome têm cobertura e granularidade diferentes; a ausência de uma via em um banco não significa ausência biológica.
- A rede é **centrada em NQO1** (primeira vizinhança). Ela não é a rede completa do interactoma humano.
- O log2FC vem do dataset GSE116250 e depende do modelo estatístico do projeto (DCM vs CTRL, ajustado por sexo/idade). Genes não medidos aparecem em cinza.
- O símbolo **SOD2-2** foi mantido como retornado pelo STRING; no cruzamento com KEGG, foi tratado como alias de **SOD2**.
- KEGG e Reactome são bases de **vias**, não de PPI. A rede PPI é do STRING.

---

## 6. Reprodutibilidade

```bash
python tutoriais/07_PPI_NQO1_string_kegg_reactome.py
```

Para forçar a atualização dos dados baixados:

```bash
NQO1_REFRESH=1 python tutoriais/07_PPI_NQO1_string_kegg_reactome.py
```

Os dados brutos ficam em `output/temp/nqo1_ppi_multidb/`. As figuras e tabelas são gravadas em `output/figures/` e `output/tables/`.

---

## 7. Arquivos gerados

### Figuras
- `output/figures/nqo1_rede_ppi_string.png`
- `output/figures/nqo1_kegg_pathways.png`
- `output/figures/nqo1_reactome_pathways.png`
- `output/figures/nqo1_heatmap_vias_parceiros.png`
- `output/figures/nqo1_kegg_diagram_hsa05208.png`
- `output/figures/nqo1_reactome_diagram_nfe2l2.png`

### Tabelas
- `output/tables/nqo1_string_ppi_arestas.csv`
- `output/tables/nqo1_string_ppi_nos.csv`
- `output/tables/nqo1_kegg_pathways.csv`
- `output/tables/nqo1_reactome_pathways.csv`
- `output/tables/nqo1_reactome_hierarquia.csv`
- `output/tables/nqo1_resumo_metricas.csv`

---

*Relatório gerado automaticamente pelo script `tutoriais/07_PPI_NQO1_string_kegg_reactome.py`.*
"""
    (ROOT / "RELATORIO_NQO1_PPI_STRING_KEGG_REACTOME.md").write_text(texto, encoding="utf-8")
    log("Relatório salvo em RELATORIO_NQO1_PPI_STRING_KEGG_REACTOME.md.")


# ---------------------------------------------------------------------------
# 9) MAIN
# ---------------------------------------------------------------------------
def main() -> None:
    log("Iniciando análise NQO1: STRING + KEGG + Reactome.")

    arestas = baixar_rede_string()
    g = montar_grafo(arestas)
    deg = carregar_deg()

    kegg, genes_kegg, _sym2id = baixar_kegg()
    reactome, hierarquia, genes_reactome, _db_id = baixar_reactome()

    parceiros = {n for n in g.nodes if n != GENE}
    kegg, reactome = completar_kegg_reactome(
        kegg, genes_kegg, reactome, genes_reactome, parceiros
    )

    nodes = anotar_nos(g, arestas, kegg, genes_kegg, reactome, genes_reactome, deg)

    # Tabelas.
    arestas.to_csv(TAB / "nqo1_string_ppi_arestas.csv", index=False)
    nodes.to_csv(TAB / "nqo1_string_ppi_nos.csv", index=False)
    kegg.to_csv(TAB / "nqo1_kegg_pathways.csv", index=False)
    reactome.to_csv(TAB / "nqo1_reactome_pathways.csv", index=False)
    hierarquia.to_csv(TAB / "nqo1_reactome_hierarquia.csv", index=False)
    log("Tabelas CSV salvas em output/tables/.")

    # Figuras.
    figura_rede(g, nodes, arestas)
    figura_kegg(kegg)
    figura_reactome(reactome)
    figura_heatmap(nodes, kegg, reactome)
    figura_kegg_diagrama()
    figura_reactome_diagrama(reactome)

    # Relatório.
    escrever_relatorio(nodes, kegg, reactome, hierarquia, arestas, g)

    log("Concluído.")
    log(f"Rede: {g.number_of_nodes()} nós, {g.number_of_edges()} arestas, "
        f"{g.degree(GENE)} parceiros diretos de NQO1.")
    log(f"KEGG: {len(kegg)} vias | Reactome: {len(reactome)} vias.")


if __name__ == "__main__":
    main()
