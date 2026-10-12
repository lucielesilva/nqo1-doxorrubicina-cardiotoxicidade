# Auditoria de Dados

> Relatório de integridade e reprodutibilidade gerado antes do versionamento.
> Projeto: `nqo1-doxorrubicina-cardiotoxicidade`
> Data da auditoria: 2026-08-12
> Ferramentas: Python 3.14 (pandas 2.3.3, numpy 2.4.6) · inspeção manual de `analise_nqo1_cardiotoxicidade.R`

---

## 0. Escopo e arquivos auditados

| Arquivo | Formato | Papel |
|---|---|---|
| `data/raw/GSE116250_rpkm.txt.gz` | TSV (gzip) | Matriz de expressão RPKM (RNA-seq) |
| `data/raw/GSE116250_series_matrix.txt.gz` | TXT (gzip) | Metadados GEO (Série) |
| `data/raw/9606.protein.{aliases,info,links}.v11.5.*.gz` | TSV (gzip) | STRING v11.5 (PPI) |
| `data/raw/GPL16791.soft.gz` | SOFT (gzip) | Anotação de plataforma |
| `output/tables/nqo1/*.csv` | CSV | Resultados NQO1 × fármacos |
| `output/tables/ampk/*.csv` | CSV | DEG / GSEA / GSVA AMPK |
| `output/tables/pubchem/*.csv` | CSV | Associação gene-composto PubChem |
| `output/tables/ppi/*.csv` | CSV | Rede PPI (nós/arestas) |
| `analise_nqo1_cardiotoxicidade.R` | R | Pipeline analítico principal |

---

## 1. Inconsistências de Tipo

### 1.1 Colunas-índice residuais (artefatos que devem ser removidos)
| Arquivo | Coluna | Problema |
|---|---|---|
| `output/tables/nqo1/dgb_resultados_nqo1_consolidados.csv` | `X1` | Índice residual **por aba** da planilha (0-1806, 1716 valores duplicados em 3523 linhas). Não carrega informação biológica; deve ser descartado na leitura. |
| `output/tables/ampk/metadados_das_amostras.csv` | `Unnamed: 0` | Cópia redundante de `geo_accession` (GSMxxxx). Remover. |

**Ação sugerida:** em `analise_nqo1_cardiotoxicidade.R`, escrever os CSVs com `row.names = FALSE` já é feito para a maioria; a coluna `X1` vem da fonte `.xlsx` (cada aba traz um índice próprio). Ao consolidar as abas, dropar a primeira coluna antes do `bind_rows`.

### 1.2 `Fold_Change` é um escore sinalizado, não uma razão
A coluna `Fold_Change` em `dgb_resultados_nqo1_consolidados.csv` **não é uma razão de fold-change clássica**. Ex.: `trichostatin A / MCF7` tem `Fold_Change = 0.459694` e, mesmo assim, é rotulada como `"Upregula NQO1"` (porque a convenção do script é `Fold_Change > 0 => up`).

> Uma razão FC = 0.46 significaria **down** (46% do controle). O uso de um valor < 1 com rótulo "up" indica que a coluna é um **escore de efeito sinalizado** (sinal = direção; magnitude = tamanho do efeito), herdado de CMAP/LINCS/CREEDS.

**Impacto:** qualquer leitor que interprete `Fold_Change` como razão concluirá o sentido errado. Recomenda-se renomear para `effect_score` (ou documentar a semântica) e, se necessário, derivar `log2FC` real por fármaco.

### 1.3 Datas em formato não-ISO
`metadados_das_amostras.csv`: `submission_date` e `last_update_date` estão como strings `"Jun 25 2018"` / `"Nov 14 2018"`. Recomenda-se `as.Date(x, format = "%b %d %Y")` para operações temporais.

### 1.4 Zeros exatos em p-value / q-value
`dgb_resultados_nqo1_consolidados.csv` contém **1 linha** com `p_value == 0` e `q_value == 0` (`trichostatin A`, aba `cmap_results_up`). p-values exatos iguais a zero são artefato de *underflow* numérico, não um valor real. O `analise_nqo1_cardiotoxicidade.R` já capa `-log10(q)` em 50 (`q_to_neglog10`), o que é correto; o valor bruto, porém, deveria ser registrado como `< 1e-300` em vez de `0` para evitar ambiguidade.

---

## 2. Incongruências de Rótulo

### 2.1 [CORRIGIDO] Fluxo analítico definitivo e divergências de rótulo
Após a auditoria forense e a clarificação do autor, o fluxo correto é:

1. **DE do transcriptoma completo** (GSE116250, RNA-seq) → NQO1 = 2º gene mais desregulado, **DOWNREGULADO** na DCM (`log2FC = −1,686`).
2. **Via AMPK (hsa04152)** avaliada por GSEA/GSVA → **não significativa** (padj = 0,419).
3. **NQO1 NÃO pertence à via AMPK**, foi selecionada a partir do transcriptoma completo, por relevância redox.
4. **DGB** (CMAP/LINCS/CREEDS) → **doxorrubicina** = downregulador de NQO1 (q = 2,21 × 10⁻²⁸).

Divergências históricas detectadas e corrigidas:
- O dataset `GSE116250` é **RNA-seq** (não "microarranjos", como aparece em versões antigas do texto).
- NQO1 é **down** (não "up") na DCM.
- NQO1 **não** foi filtrada a partir da via AMPK.
- Não há rótulo `Doxo` × `Controle` no dataset de expressão; a associação doxorrubicina↔NQO1 vem das bases DGB (LINCS L1000 + CREEDS).

**Conclusão:** a hipótese de cardiotoxicidade é **indireta/integrativa** e deve ser assim apresentada.

### 2.2 Deriva de rótulo do grupo controle (3 nomes para o mesmo grupo)
| Fonte | Rótulo do grupo controle |
|---|---|
| `metadados_das_amostras.csv` (`disease:ch1`) | `non-failing` |
| Matriz RPKM (prefixo de coluna) | `NF` |
| `expressao_diferencial_completa.csv` (contraste) | controle implícito = NF |
| `gsva_via_ampk.csv` (`group`) | `CTRL` |

Um mesmo grupo é chamado de `non-failing`, `NF` e `CTRL`. Recomenda-se padronizar (ex.: `CTRL` em todos os artefatos derivados) e registrar a correspondência no README.

### 2.3 Escopo divergente entre tabelas
- `expressao_diferencial_completa.csv` cobre as 64 amostras (NF + DCM + ICM).
- `gsva_via_ampk.csv` contém **51 amostras** (37 DCM + 14 CTRL), **excluindo ICM**. A exclusão de ICM é plausível (foco DCM), mas deve ser declarada; hoje não está documentada no script.

### 2.4 Consistência de identificadores (positivo)
A conferência direta entre os **64 títulos** de `metadados_das_amostras.csv` e as **64 colunas** de amostra da matriz RPKM resultou em **conjuntos idênticos** (0 amostras órfãs, 0 ausentes). ✔️ Não há amostra com metadado sem expressão nem o inverso.

---

## 3. Dados Faltantes (Missing)

| Conjunto | Variável | % NA | Alerta (>5%)? |
|---|---|---|---|
| RPKM (GSE116250) | todas as células (3.710.336) | **0,0000%** | Não |
| `expressao_diferencial_completa.csv` | todas | **0%** | Não |
| `expressao_diferencial_via_ampk.csv` | todas | **0%** | Não |
| `dgb_resultados_nqo1_consolidados.csv` | `CREEDS_ID` | 96,48% | ⚠️ informativo* |
| `dgb_resultados_nqo1_consolidados.csv` | `GEO_ID` | 96,48% | ⚠️ informativo* |
| `dgb_resultados_nqo1_consolidados.csv` | `DrugBank_ID` | 98,18% | ⚠️ informativo* |
| `dgb_resultados_nqo1_consolidados.csv` | `LINCS_sig_id`/`LINCS_pert_id` | 4,09% | Não |
| `dgb_resultados_nqo1_consolidados.csv` | `Cell_Line`/`Time`/`Dose` | 3,52% | Não |
| `pubchem_associacoes_compostos.csv` | `Name` | **100%** | ⚠️ esperado** |

\* As colunas de ID específico de banco (`CREEDS_ID`, `GEO_ID`, `DrugBank_ID`) são naturalmente vazias para registros oriundos de outros bancos. Não é erro, mas confirma que o `NA` está sendo usado como "não se aplica". Recomenda-se distinguir `NA` (ausente) de "não se aplica" (ex.: string `"n/a"`), se isso importar para análises futuras.

\** `pubchem_associacoes_compostos.csv` tem `Name` 100% NA porque o nome só é resolvido na versão `pubchem_associacoes_compostos_com_nomes.csv`. Recomenda-se fundir num único artefato para não manter uma coluna permanentemente vazia.

**Conclusão de missing:** nenhuma variável analítica contínua (expressão, FC, p/q-value) ultrapassa 5% de NA. ✔️

---

## 4. Duplicatas

| Tabela | Linhas | Duplicadas (linhas) | Duplicadas (chave) |
|---|---|---|---|
| `dgb_resultados_nqo1_consolidados.csv` | 3.523 | 0 | `X1` tem 1.716 duplicados (artefato, ver §1.1) |
| `expressao_diferencial_completa.csv` | 41.842 | 0 | `gene` único (0 duplicados) |
| `expressao_diferencial_via_ampk.csv` | 120 | 0 | `gene` único |
| `metadados_das_amostras.csv` | 64 | 0 | `geo_accession`/`title` únicos |
| `ppi_ampk_string_arestas.csv` | 1.874 | 0 |, |
| `ppi_ampk_string_nos.csv` | 119 | 0 | `gene` único |
| `gsea_todas_as_vias_kegg.csv` | 187 | 0 | `pathway` único |

**Conclusão:** nenhuma linha duplicada nas tabelas primárias. ✔️

---

## 5. Benchmarks e Thresholds (Validação Estatística/Biológica)

### 5.1 Threshold de Viabilidade Celular (MTT/SRB)
- **Status: N/A, não há ensaio de viabilidade no projeto.**
- Nenhum CSV contém absorbância/viabilidade normalizada. O efeito da doxorrubicina é avaliado **indiretamente** pela assinatura de expressão de NQO1 (`mean_fc = -1.065`, `min_q = 2.21e-28`, 15 testes em LINCS L1000 + CREEDS).
- **Observação:** a validação in vitro está fora do escopo deste projeto, que é restrito a análises in silico com dados públicos. Não há threshold de viabilidade a avaliar.

### 5.2 Threshold de Expressão Gênica (NQO1)
- **NQO1 no DEG (GSE116250, DCM vs controle):** `log2FoldChange = -1.6856`, `padj = 6.86e-13` → **PASSA** no benchmark. Em escala linear, `FC = 2^-1.6856 ≈ 0.31` (≈ 3,2× de redução), muito além do limiar `|log2FC| = 0.585` (= log2 1.5). ✔️
- **Doxorrubicina → NQO1 (assinatura de fármaco):** `mean_fc = -1.065`, `median_fc = -0.807`, `max_abs_fc = 4.21` → **AVISO**. O efeito é estatisticamente significativo (q = 2.21e-28), mas o **tamanho de efeito médio está abaixo de 1,5×** (|−1.065| < 1.5). A significância é real, porém o efeito médio é heterogêneo (doses/linhagens distintas). Reportar como "significativo, porém com efeito médio < 1,5×, heterogêneo entre condições".

### 5.3 Outliers (Z-score, |Z| > 3.0)
| Variável | n outliers (|Z|>3) |
|---|---|
| `Fold_Change` | 6 |
| `Specificity` | 29 |
| `p_value` | 99 |
| `q_value` | 4 |
| `score` | 78 |
| `neg_log10_q_capped` | 77 |

- Os outliers em `p_value`/`q_value`/`score`/`neg_log10_q_capped` são **esperados** em dados de assinatura em larga escala (q-values extremamente pequenos). Não indicam erro de coleta.
- **Ação:** aplicar winsorização/estatística robusta na variável `score` (ranking de fármacos) e preservar os valores brutos em `/data/raw` (backup). Não remover linhas sem justificativa biológica.

### 5.4 Benchmark de Normalização
- A matriz `GSE116250_rpkm.txt.gz` está em **RPKM linear** (Cufflinks), **não** log2/TMM.
- O `expressao_diferencial_completa.csv` foi gerado por limma (colunas `AveExpr`, `B`, `t`), o que pressupõe dados em escala **log2**, porém a **proveniência da normalização não está documentada** (não há registro de `voom`/`TMM`/`CPM` no repositório).
- **Ação:** registrar no README/script que o DEG foi computado sobre `log2(CPM)` via `edgeR::voom` (ou método equivalente), e adicionar `log2(RPKM + 1)` como transformação reprodutível antes do DE. O `config.yaml` define `normalization.expected: log2-TMM`.

### 5.5 Resultado-chave a não superestimar: GSEA AMPK
- `gsea_via_ampk.csv`: **`NES = -1.193`, `padj = 0.419`** → **NÃO significativo** (FDR > 0.05).
- O enriquecimento da via AMPK **não sobrevive** à correção de múltiplas comparações. O README e qualquer manuscrito derivado devem evitar afirmar "ativação/inibição significativa da via AMPK" a partir deste teste.

---

## 6. Sumário Executivo da Auditoria

| Item | Veredito |
|---|---|
| Missing em dados analíticos | ✅ OK (< 5%) |
| Duplicatas em tabelas primárias | ✅ OK (0) |
| Consistência amostra × metadado | ✅ OK (64 = 64, conjuntos idênticos) |
| Colunas-índice residuais (`X1`, `Unnamed: 0`) | ⚠️ remover |
| Semântica de `Fold_Change` (escore vs razão) | ⚠️ documentar/renomear |
| Zeros exatos de p/q-value | ⚠️ tratar como underflow |
| Rótulo do grupo controle (3 nomes) | ⚠️ padronizar |
| Escopo GSVA (ICM excluído) | ⚠️ documentar |
| Narrativa "doxo" × dado real (DCM/ICM) | 🔴 CRÍTICO, esclarecer no README |
| NQO1 DEG (log2FC −1.69, padj 6.9e-13) | ✅ passa no FC>1.5 |
| Doxo→NQO1 (mean_fc −1.07, q 2.2e-28) | 🟡 AVISO (efeito médio < 1.5×) |
| Viabilidade MTT/SRB | ⬜ N/A (dado ausente) |
| GSEA AMPK (NES −1.19, padj 0.42) | 🔴 não significativo |
| Normalização documentada | ⚠️ adicionar proveniência |

**Recomendação final:** os dados estão estruturalmente íntegros e sem duplicatas/missing críticos; os problemas são majoritariamente de **documentação e rótulo**, não de corrupção. Antes da publicação, resolver os itens 🔴/⚠️ acima (narrativa, renomear `Fold_Change`, padronizar `CTRL`, documentar normalização).
