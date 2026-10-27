# NQO1 — rede PPI no STRING e vias no KEGG e Reactome

**Gene alvo:** NQO1 (UniProt P15559) · **Espécie:** *Homo sapiens* (9606)  
**Dataset de expressão:** GSE116250 (DCM vs CTRL, limma)  
**Data da análise:** 26/10/2026

---

## 1. Resumo

| Banco | O que foi consultado | Resultado principal |
|---|---|---|
| **STRING** | Rede de vizinhança direta de NQO1, score ≥ 700 | **43 nós** e **221 arestas**; NQO1 tem **42 parceiros diretos** |
| **KEGG** | Vias que contêm NQO1 e genes dessas vias | **7 vias**; NQO1 aparece em vias de metabolismo e de doenças humanas |
| **Reactome** | Vias de baixo nível que contêm NQO1 e hierarquia | **3 vias**; a principal é a via **NFE2L2/KEAP1** |

> **NQO1 no dataset:** log2FC = **-1.66** (padj = **6.50e-13**), ou seja, **downregulado** na cardiomiopatia dilatada (DCM) em relação ao controle.

![Rede PPI centrada em NQO1](output/figures/nqo1_rede_ppi_string.png)

**Como ler a rede:** cada bolinha é uma proteína; o centro é o NQO1. O anel externo mostra os 42 parceiros diretos. O tamanho do nó representa o grau (nº de conexões na rede); a cor representa o log2FC na DCM (azul = down, vermelho = up; cinza = gene não medido no dataset). As linhas escuras ligam NQO1 aos parceiros; as linhas claras mostram interações entre os próprios parceiros.

---

## 2. Métodos (resumo)

1. **STRING** — API `network` (`species=9606`, `required_score=700`, `limit=500`). O score combinado varia de 0 a 1; ≥ 0,700 é considerado alta confiança. A rede inclui NQO1, seus parceiros diretos e as interações entre esses parceiros.
2. **KEGG** — API REST (`/link/pathway/hsa:1728`, `/link/hsa/{pathway}`, `/get/{pathway}`). Foram obtidas as vias que contêm NQO1, seus nomes, classes e genes.
3. **Reactome** — Content Service (`/search/query`, `/data/pathways/low/entity/...`, `/data/event/.../ancestors`, `/data/participants/...`). Foram obtidas as vias de baixo nível, a hierarquia até o top-level e os participantes.
4. **Expressão** — resultado do `limma` do projeto (`resultado_limma_DCM_vs_CTRL.rds`), usado apenas para colorir os nós.
5. **Cruzamento** — cada parceiro do STRING foi comparado com os genes/participantes das vias KEGG e Reactome que contêm NQO1.

> **Importante:** KEGG e Reactome **não são bases de PPI**. A rede de interações proteína-proteína vem do STRING; KEGG/Reactome entram como anotação de vias e contexto funcional.

---

## 3. Resultados

### 3.1 Rede PPI do STRING

A rede tem **43 proteínas** e **221 interações** de alta confiança. O NQO1 é o nó central, com **42 parceiros diretos**. Os nós mais conectados (hubs) são **NFE2L2, KEAP1, HMOX1, TXN, GSR**.

**Tabela 1 — 15 parceiros com maior grau na rede** (o NQO1 não entra nesta lista; ele é o centro).

| Gene   |   Grau |   Score NQO1 | log2FC   | padj     |   Vias KEGG* |   Vias Reactome* |
|:-------|-------:|-------------:|:---------|:---------|-------------:|-----------------:|
| NFE2L2 |     25 |        0.984 | +0.11    | 2.83e-01 |            4 |                2 |
| KEAP1  |     21 |        0.958 | -0.26    | 3.41e-04 |            4 |                0 |
| HMOX1  |     20 |        0.915 | -0.00    | 9.96e-01 |            5 |                1 |
| TXN    |     19 |        0.702 | -0.38    | 8.46e-04 |            1 |                1 |
| GSR    |     17 |        0.816 | -0.29    | 1.16e-03 |            1 |                1 |
| GPX2   |     17 |        0.772 | —        | —        |            1 |                0 |
| SOD2-2 |     17 |        0.786 | —        | —        |            1 |                0 |
| GPX8   |     16 |        0.735 | +0.52    | 3.77e-02 |            1 |                0 |
| GPX5   |     16 |        0.719 | —        | —        |            1 |                0 |
| GPX6   |     16 |        0.719 | —        | —        |            1 |                0 |
| GPX7   |     16 |        0.72  | +0.43    | 6.22e-04 |            1 |                0 |
| GPX3   |     16 |        0.727 | -0.67    | 3.67e-05 |            1 |                0 |
| GCLC   |     14 |        0.928 | -0.37    | 4.05e-04 |            2 |                1 |
| TNF    |     14 |        0.799 | —        | —        |            1 |                0 |
| TXNRD1 |     13 |        0.774 | -0.81    | 4.51e-04 |            2 |                1 |

\* Número de vias **de NQO1** (KEGG ou Reactome) em que o parceiro também aparece; não é o total de vias do gene.

**Arquivos completos:**
- `output/tables/nqo1_string_ppi_arestas.csv` — todas as arestas, com score e canais de evidência do STRING (pode ser importado no Cytoscape como edge list).
- `output/tables/nqo1_string_ppi_nos.csv` — todos os nós, com grau, centralidades, log2FC e vias compartilhadas.

#### Módulos que aparecem na rede

- **Módulo NRF2/KEAP1:** NFE2L2, KEAP1, GCLC, GCLM, HMOX1, GSR, TXNRD1, SRXN1, TXN. É o eixo clássico de resposta antioxidante e detoxificação.
- **Família GPX:** GPX2, GPX3, GPX5, GPX6, GPX7 e GPX8 formam um grupo muito conectado, ligado ao metabolismo de peróxidos.
- **Eixo TP53/STUB1/PARK7:** TP53, TP73, STUB1 e PARK7 aparecem como conexões de estresse celular e controle de qualidade proteica.
- **Metabolismo de quinonas/vitamina K:** GGCX, VKORC1, VKORC1L1, CRYZ e NDUFS7.

### 3.2 KEGG

O NQO1 está anotado em **7 vias KEGG**:

![Vias KEGG que contêm NQO1](output/figures/nqo1_kegg_pathways.png)

**Tabela 2 — Vias KEGG que contêm NQO1.**

| ID       | Via KEGG                                            | Classe                                           |   Genes na via |   Parceiros STRING | Parceiros                                                                                                           |
|:---------|:----------------------------------------------------|:-------------------------------------------------|---------------:|-------------------:|:--------------------------------------------------------------------------------------------------------------------|
| hsa01100 | Metabolic pathways                                  | Metabolism; Global and overview maps             |           1583 |                 21 | ADK;CYP1A1;GCLC;GCLM;GGCX;GOT2;GPX2;GPX3;GPX5;GPX6;GPX7;GPX8;GSR;GSTM1;GSTP1;HMOX1;NDUFS7;NME1;ODC1;VKORC1;VKORC1L1 |
| hsa05208 | Chemical carcinogenesis - reactive oxygen species   | Human Diseases; Cancer: overview                 |            226 |                 11 | AKT1;CYP1A1;EPHX1;GSTM1;HIF1A;HMOX1;JUN;KEAP1;NDUFS7;NFE2L2;SOD2-2                                                  |
| hsa05418 | Fluid shear stress and atherosclerosis              | Human Diseases; Cardiovascular disease           |            142 |                 10 | AKT1;GSTM1;GSTP1;HMOX1;JUN;KEAP1;NFE2L2;TNF;TP53;TXN                                                                |
| hsa05200 | Pathways in cancer                                  | Human Diseases; Cancer: overview                 |            533 |                 10 | AKT1;GSTM1;GSTP1;HIF1A;HMOX1;JUN;KEAP1;NFE2L2;TP53;TXNRD1                                                           |
| hsa05225 | Hepatocellular carcinoma                            | Human Diseases; Cancer: specific types           |            170 |                  8 | AKT1;GSTM1;GSTP1;HMOX1;KEAP1;NFE2L2;TP53;TXNRD1                                                                     |
| hsa01240 | Biosynthesis of cofactors                           | Metabolism; Metabolism of cofactors and vitamins |            154 |                  6 | GCLC;GCLM;GGCX;NME1;VKORC1;VKORC1L1                                                                                 |
| hsa00130 | Ubiquinone and other terpenoid-quinone biosynthesis | Metabolism; Metabolism of cofactors and vitamins |             12 |                  3 | GGCX;VKORC1;VKORC1L1                                                                                                |

**Como ler:** a barra mostra quantos parceiros diretos de NQO1 no STRING também estão naquela via. O número entre parênteses é o total de genes da via no KEGG.

**Destaque:** `hsa05208` (*Chemical carcinogenesis — reactive oxygen species*) é particularmente relevante para o projeto, pois conecta NQO1 ao estresse oxidativo e à sinalização NRF2/KEAP1. O mapa dessa via, com o NQO1 destacado em vermelho, está em:

<img src="output/figures/nqo1_kegg_diagram_hsa05208.png" width="720" alt="Mapa KEGG hsa05208 com NQO1 destacado">

*No mapa, os retângulos vermelhos marcam as ocorrências de NQO1.*

### 3.3 Reactome

O NQO1 está anotado em **3 vias de baixo nível** no Reactome:

![Vias Reactome que contêm NQO1](output/figures/nqo1_reactome_pathways.png)

**Tabela 3 — Vias Reactome de baixo nível que contêm NQO1.**

| ID            | Via Reactome                                          | Top-level                     |   Genes/participantes |   Parceiros STRING | Parceiros                                        |
|:--------------|:------------------------------------------------------|:------------------------------|----------------------:|-------------------:|:-------------------------------------------------|
| R-HSA-9818027 | NFE2L2 regulating anti-oxidant/detoxification enzymes | Cellular responses to stimuli |                    20 |                  9 | GCLC;GCLM;GSR;HMOX1;MAFK;NFE2L2;SRXN1;TXN;TXNRD1 |
| R-HSA-9943962 | CHD6, CHD7, CHD8, CHD9 subfamily                      | Chromatin organization        |                    49 |                  2 | MAFK;NFE2L2                                      |
| R-HSA-350562  | Regulation of ornithine decarboxylase (ODC)           | Metabolism                    |                    38 |                  1 | ODC1                                             |

**Hierarquia até o top-level** (eventos ancestrais, sem repetir a própria via):

| Via NQO1      | Evento ancestral                          | Top-level                     |
|:--------------|:------------------------------------------|:------------------------------|
| R-HSA-350562  | Metabolism of polyamines                  | Metabolism                    |
| R-HSA-350562  | Metabolism of amino acids and derivatives | Metabolism                    |
| R-HSA-350562  | Metabolism                                | Metabolism                    |
| R-HSA-9943962 | CHD chromatin remodelers                  | Chromatin organization        |
| R-HSA-9943962 | ATP-dependent chromatin remodelers        | Chromatin organization        |
| R-HSA-9943962 | Chromatin organization                    | Chromatin organization        |
| R-HSA-9818027 | Nuclear events mediated by NFE2L2         | Cellular responses to stimuli |
| R-HSA-9818027 | KEAP1-NFE2L2 pathway                      | Cellular responses to stimuli |
| R-HSA-9818027 | Cellular response to chemical stress      | Cellular responses to stimuli |
| R-HSA-9818027 | Cellular responses to stress              | Cellular responses to stimuli |
| R-HSA-9818027 | Cellular responses to stimuli             | Cellular responses to stimuli |

**Destaque:** a via **R-HSA-9818027** (*NFE2L2 regulating anti-oxidant/detoxification enzymes*) é a que mais converge com a rede PPI: **9 parceiros** do STRING também participam dela. Ela sobe na hierarquia até **Cellular responses to stimuli**, passando por *KEAP1-NFE2L2 pathway* e *Cellular response to chemical stress*.

<img src="output/figures/nqo1_reactome_diagram_nfe2l2.png" width="720" alt="Diagrama Reactome da via NFE2L2">

*Diagrama da via R-HSA-9818027 no Reactome.*

### 3.4 Integração: parceiros que compartilham vias com NQO1

![Heatmap de vias compartilhadas](output/figures/nqo1_heatmap_vias_parceiros.png)

O heatmap mostra, para cada parceiro de NQO1, em quais vias KEGG (bloco à esquerda) e Reactome (bloco à direita) ele também aparece. As linhas mais “cheias” são os parceiros que mais convergem funcionalmente com NQO1. O NQO1 aparece em todas as colunas por definição (é o gene consultado).

**Tabela 4 — Parceiros com maior número de vias compartilhadas.**

| Gene     |   Grau |   KEGG |   Reactome |   Total | IDs KEGG                                     | IDs Reactome                |
|:---------|-------:|-------:|-----------:|--------:|:---------------------------------------------|:----------------------------|
| NFE2L2   |     25 |      4 |          2 |       6 | hsa05200;hsa05208;hsa05225;hsa05418          | R-HSA-9818027;R-HSA-9943962 |
| HMOX1    |     20 |      5 |          1 |       6 | hsa01100;hsa05200;hsa05208;hsa05225;hsa05418 | R-HSA-9818027               |
| GSTM1    |     11 |      5 |          0 |       5 | hsa01100;hsa05200;hsa05208;hsa05225;hsa05418 |                             |
| KEAP1    |     21 |      4 |          0 |       4 | hsa05200;hsa05208;hsa05225;hsa05418          |                             |
| GSTP1    |     13 |      4 |          0 |       4 | hsa01100;hsa05200;hsa05225;hsa05418          |                             |
| AKT1     |     12 |      4 |          0 |       4 | hsa05200;hsa05208;hsa05225;hsa05418          |                             |
| GCLC     |     14 |      2 |          1 |       3 | hsa01100;hsa01240                            | R-HSA-9818027               |
| TXNRD1   |     13 |      2 |          1 |       3 | hsa05200;hsa05225                            | R-HSA-9818027               |
| TP53     |     12 |      3 |          0 |       3 | hsa05200;hsa05225;hsa05418                   |                             |
| GCLM     |      9 |      2 |          1 |       3 | hsa01100;hsa01240                            | R-HSA-9818027               |
| JUN      |      7 |      3 |          0 |       3 | hsa05200;hsa05208;hsa05418                   |                             |
| GGCX     |      4 |      3 |          0 |       3 | hsa00130;hsa01100;hsa01240                   |                             |
| VKORC1   |      4 |      3 |          0 |       3 | hsa00130;hsa01100;hsa01240                   |                             |
| VKORC1L1 |      3 |      3 |          0 |       3 | hsa00130;hsa01100;hsa01240                   |                             |
| TXN      |     19 |      1 |          1 |       2 | hsa05418                                     | R-HSA-9818027               |

---

## 4. Métricas do resumo

| Métrica                              | Valor               |
|:-------------------------------------|:--------------------|
| Gene alvo                            | NQO1                |
| UniProt                              | P15559              |
| Espécie                              | Homo sapiens (9606) |
| Rede STRING — nós                    | 43                  |
| Rede STRING — arestas                | 221                 |
| Parceiros diretos de NQO1            | 42                  |
| Score mínimo STRING                  | 700                 |
| Vias KEGG com NQO1                   | 7                   |
| Vias Reactome (baixo nível) com NQO1 | 3                   |
| NQO1 log2FC (DCM vs CTRL)            | -1.664              |
| NQO1 padj                            | 6.50e-13            |
| NQO1 direção                         | down                |

---

## 5. Limitações

- O STRING integra **evidências heterogêneas** (texto, experimentos, bancos de dados, coexpressão, vizinhança genômica, fusão gênica). Uma aresta não significa necessariamente contato físico direto.
- O score ≥ 700 reduz falsos positivos, mas ainda pode incluir associações funcionais indiretas.
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
