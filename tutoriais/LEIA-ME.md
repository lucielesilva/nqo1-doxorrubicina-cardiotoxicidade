# NQO1 — roteiro de aula (scripts separados)

Esta pasta reúne a aula em **etapas numeradas**. O pipeline original
(`analise_nqo1_cardiotoxicidade.R`, na raiz do projeto) foi quebrado em scripts
curtos, cada um com muitos comentários explicando o que está sendo feito e por
quê.

A regra da pasta é simples: **o número do arquivo é a ordem de execução**.

## Conteúdo da pasta

| Arquivo | Tipo |
|---|---|
| `00_configuracao_e_dados.R` … `06_batch_deconvolution.R` | scripts da aula (R) |
| `APRESENTACAO_NQO1_TUTORIAIS.pptx` | apresentação da aula (31 slides, editável) |
| `montar_apresentacao.py` | script que gera a apresentação a partir dos resultados |
| `LEIA-ME.md` | este arquivo |

## Apresentação

`APRESENTACAO_NQO1_TUTORIAIS.pptx` é a apresentação usada em sala: **50 slides**,
formato 16:9, com **texto, tabelas e figuras editáveis**, texto alternativo nas
figuras e notas do apresentador em cada slide. Ela abre normalmente no
PowerPoint, no LibreOffice Impress e no **Google Apresentações** (basta importar
o arquivo no Drive).

Blocos da apresentação: objetivos e glossário → contexto clínico (DOX,
mecanismos, morte celular, NQO1, cardioproteção) → dados e métodos → resultados
transcriptômicos (UMAP, volcano, NQO1, heatmap) → do gene ao sistema → via e
robustez → farmacogenômica → conclusões, discussão e referências.

Os slides de mecanismos usam as figuras dos artigos em `tutoriais/artigos/`:
Kong et al. (2022), Fig. 1 (ciclo redox e estresse oxidativo) e Fig. 2
(apoptose, piroptose e ferroptose); Rodrigues De Souza et al. (2025), Fig. 1
(as 25 formas de morte celular regulada) e Fig. 3 (modelo de multioptose).

Os números e as figuras dos slides vêm das tabelas em `output/tables/` e das
imagens em `output/figures/`. Se você rodar a análise de novo e os valores
mudarem, gere a apresentação outra vez com:

```powershell
python tutoriais\montar_apresentacao.py
```

O script lê os CSVs de resultado, recalcula os números e monta os slides (título,
tópicos, tabelas, figuras, notas e numeração).

## Ordem dos scripts

| Script | O que faz | Precisa de qual antes? |
|---|---|---|
| `00_configuracao_e_dados.R` | Descobre a raiz do projeto, cria as pastas de saída, baixa o GEO e a via AMPK, e guarda um cache. **Não faz análise.** | — |
| `01_UMAP.R` | UMAP das amostras (exploração não supervisionada). | 00 |
| `02_volcano.R` | limma DCM vs CTRL (ajustado por sexo/idade) + tabela de expressão + **volcano plot**. **A aula pode terminar aqui.** | 00 |
| `03_gene_sets.R` | Separa os genes do dataset dos genes da via AMPK (Venn e tabela). | 02 (gera o resultado do limma) |
| `04_rede_PPI.R` | Rede de interações da via AMPK no STRING. | 02 |
| `05_extras_GSEA_GSVA.R` | GSEA e GSVA da via AMPK. | 02 |
| `06_batch_deconvolution.R` | Correção de batch com SVA e deconvolução celular com xCell. | 02 |
| `executar_tutoriais.R` (na raiz) | Runner que executa tudo na ordem correta. | — |

O `03_gene_sets.R` vem antes da rede e do GSEA de propósito: primeiro fica claro
*quais* genes entraram na análise, depois se interpretam vias e interações.

## Como rodar

### Pipeline completo

Da raiz do projeto:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' 'executar_tutoriais.R'
```

O runner executa, nesta ordem: `01_UMAP` → `02_volcano` → `03_gene_sets` →
`04_rede_PPI` → `05_extras_GSEA_GSVA` → `06_batch_deconvolution`, carregando o
`00` uma única vez.

### Aula “até o volcano”

```powershell
cd "C:\Users\oorie\OneDrive\Documentos\TRABALHOS\nqo1-doxorrubicina-cardiotoxicidade"
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' 'tutoriais\01_UMAP.R'
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' 'tutoriais\02_volcano.R'
```

### Um script isolado

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' 'tutoriais\02_volcano.R'
```

Cada script de análise carrega o `00` sozinho (`source()`), então não é preciso
rodar o `00` antes — a não ser que você queira baixar os dados antecipadamente.

## Onde os arquivos ficam

Os scripts estão em `tutoriais/`, mas **os dados e os resultados ficam na raiz do
projeto**:

```
data/raw/        arquivos baixados do GEO e cache do STRING
data/processed/  objetos .rds (dados tratados e resultado do limma)
output/tables/   tabelas .csv
output/figures/  figuras .png
```

O `00` descobre a raiz subindo um nível a partir da própria pasta (procura o
`DESCRIPTION` do projeto, ou as pastas `data/` e `output/`). Assim, rodar um
script de dentro de `tutoriais/` ou da raiz grava sempre nos mesmos lugares —
não aparecem cópias de `data/` e `output/` dentro de `tutoriais/`.

## Por que existe o script 00?

Porque `01`, `02`, `03`, `04`, `05` e `06` precisam dos **mesmos dados**. O `00`
baixa uma vez e salva um cache em:

```
data/processed/dados_GSE116250.rds
```

Assim a aula não fica lenta baixando o GEO várias vezes. Para forçar um novo
download, apague esse arquivo ou chame `preparar_dados(forcar = TRUE)`.

## O que cada script gera

### `00_configuracao_e_dados.R`
- `output/tables/metadados_das_amostras.csv`
- `output/tables/genes_da_via_ampk_kegg.csv`
- cache dos dados em `data/processed/`

### `01_UMAP.R`
- `output/figures/UMAP_amostras.png`
- `output/tables/UMAP_coordenadas.csv`
- `output/tables/UMAP_resumo_grupos.csv`

### `02_volcano.R`
- `output/figures/volcano_plot.png`
- `output/tables/expressao_diferencial_completa.csv`
- `output/tables/expressao_nqo1_diferencial.csv`
- `data/processed/resultado_limma_DCM_vs_CTRL.rds` (usado pelos scripts 03–06)

### `03_gene_sets.R`
- `output/figures/Venn_dataset_AMPK.png`
- `output/tables/genes_dataset_vs_ampk.csv`

### `04_rede_PPI.R`
- `output/figures/rede_ppi_ampk.png`
- `output/tables/expressao_diferencial_via_ampk.csv`
- `output/tables/ppi_ampk_string_arestas.csv`
- `output/tables/ppi_ampk_string_nos.csv`
- `output/tables/ppi_ampk_hubs.csv`

A rede usa layout force-directed com mais iterações, linha mais fina/transparente,
espessura da aresta proporcional à confiança e rótulos com `ggrepel` apenas nos
20 hubs (top 5 em negrito). Isso reduz o efeito de "novelo" da rede.

### `05_extras_GSEA_GSVA.R`
- `output/figures/GSEA_AMPK_enrichment_plot.png`
- `output/figures/GSVA_AMPK_boxplot.png`
- `output/tables/gsea_via_ampk.csv`
- `output/tables/gsea_via_ampk_leading_edge.csv`
- `output/tables/gsva_via_ampk.csv`
- `output/tables/gsva_via_ampk_resumo_grupos.csv`
- `output/tables/gsva_via_ampk_teste.csv` (p-valor e `-log10(p)`)

O script está dividido em etapas explícitas: montar o ranking, rodar o GSEA,
organizar a tabela, desenhar o gráfico; depois preparar a matriz, rodar a GSVA,
resumir por grupo, testar com Wilcoxon e desenhar o boxplot.

### `06_batch_deconvolution.R`
- `output/tables/SVA_comparacao.csv`
- `output/tables/expressao_diferencial_com_SVA.csv`
- `output/tables/SVA_fatores.csv`
- `output/figures/SVA_logFC_antes_depois.png`
- `output/tables/xCell_scores.csv`
- `output/tables/xCell_teste.csv`
- `output/figures/xCell_heatmap.png`
- `output/figures/xCell_boxplots.png`

## Dependências

Pacotes do CRAN:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "readr", "igraph",
                     "ggrepel", "uwot", "ggforce", "pheatmap", "rstudioapi"))
```

Pacotes do Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "GEOquery", "Biobase", "limma", "KEGGREST",
  "fgsea", "GSVA", "STRINGdb", "sva"
))

# xCell é instalado a partir do GitHub:
# remotes::install_github("dviraran/xCell")
```

O `STRINGdb` é opcional: se ele não estiver instalado, o `04_rede_PPI.R` tenta
usar a API pública do STRING.

## Resultado rápido (dados já processados neste projeto)

- **NQO1:** log2FC = **-1,664** | padj = **6,5e-13** → downregulado na DCM.
- **UMAP:** o grupo CTRL forma um agrupamento separado; DCM e ICM se sobrepõem. O resultado depende de `n_neighbors`, `min_dist` e do número de genes.
- **Rede PPI:** **105 nós** e **837 arestas** (score >= 700).
  Top hubs: AKT1, MTOR, PRKAA1, PRKAA2, AKT2.
- **GSEA AMPK:** NES = **-1,062** | p = **0,307** (sem enriquecimento
  significativo neste recorte).
- **GSVA AMPK:** Wilcoxon DCM vs CTRL **p = 0,000287** | **-log10(p) = 3,54**.
- **SVA:** NQO1 permanece downregulado (log2FC ≈ -1,65; padj ≈ 7,4e-13).
- **xCell:** miócitos/endotélio linfático ↑; macrófagos M2/M1 ↓ (escores relativos).

Esses números podem mudar se você alterar filtros, limiares ou a versão dos
pacotes. Eles servem como referência para a aula.

## Dicas para a aula

1. Abra os scripts na ordem numérica. Cada um começa com uma pergunta e termina
   com "FIM DO SCRIPT".
2. No `02_volcano.R`, pare no volcano. As seções seguintes (gene sets, PPI,
   GSEA/GSVA) são a continuação da história, não pré-requisitos.
3. Os comentários usam analogias (fotografia, mapa, rede social) para explicar
   UMAP, volcano e PPI.
4. Se a internet falhar no meio, os caches de `data/` e `output/tables/`
   permitem continuar sem baixar tudo de novo.

## Histórico da organização

- O `01_PCA.R` foi **removido**: o UMAP (`01_UMAP.R`) é a exploração padrão do
  roteiro e o PCA não era mais usado na aula.
- A numeração foi realinhada à ordem de execução. O antigo `06_gene_sets.R`
  virou `03_gene_sets.R`, `03_rede_PPI.R` virou `04_rede_PPI.R`,
  `04_extras_GSEA_GSVA.R` virou `05_extras_GSEA_GSVA.R` e
  `05_batch_deconvolution.R` virou `06_batch_deconvolution.R`.
- A detecção da raiz do projeto passou a subir um nível a partir de `tutoriais/`,
  para que `data/` e `output/` fiquem sempre na raiz.
