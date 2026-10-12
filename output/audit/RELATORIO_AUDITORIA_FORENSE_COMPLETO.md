# 🔬 RELATÓRIO DE AUDITORIA CIENTÍFICA FORENSE

**Projeto**: resumo_expandido_luciele_cardiologia_2026
**Data da auditoria**: 2026-06-26
**Auditor**: Revisão forense completa, 15 etapas
**Arquivos analisados**: script.R, resumo_expandido .docx, todos os outputs (.csv, .xlsx, .png, .md, .RData), diretórios data/ e output/

---

## ETAPA 1, RECONSTRUÇÃO DO FLUXO CIENTÍFICO REAL

### Fluxograma do pipeline REAL (o que o script.R efetivamente faz)

```
┌─────────────────────────────────────────────────────────────────────┐
│                     FLUXO REAL DO script.R                           │
├─────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  1. INICIALIZAÇÃO                                                    │
│     ├── Define diretórios (data/, output/)                          │
│     ├── Cria estrutura de pastas                                     │
│     ├── Inicializa AUDIT_EVENTS (data.frame vazio)                   │
│     └── Carrega pacotes: dplyr, ggplot2, stringr, openxlsx, igraph  │
│                                                                      │
│  2. LIMPEZA DE ARQUIVOS OBSOLETOS                                    │
│     └── Remove arquivos antigos (PCA, boxplot legado, etc.)         │
│                                                                      │
│  3. LEITURA DO XLSX PRONTO (DGB_results_NQO1.xlsx)                  │
│     ├── 6 abas: cmap/l1000/creeds × up/down                         │
│     ├── Padroniza nomes de colunas                                   │
│     ├── Converte decimais (vírgula → ponto)                         │
│     ├── Calcula score = -log10(q) × |Fold_Change|                   │
│     ├── Classifica direção (up/down/sem direção)                    │
│     └── Filtra 3523 linhas válidas, 1281 fármacos únicos            │
│                                                                      │
│  4. RANKING DE FÁRMACOS                                              │
│     ├── Agrupa por Drug_Name (CASE SENSITIVE ⚠️)                     │
│     ├── Ordena por min_q_value, score, max_abs_fc                    │
│     └── Top 20 salvo                                                   │
│                                                                      │
│  5. FIGURAS NQO1                                                     │
│     ├── Volcano plot (Fold Change × -log10(q))                      │
│     └── Dose-response (Fold Change × dose µM)                       │
│                                                                      │
│  6. PUBCHEM (PRONTO)                                                 │
│     ├── Carrega pubchem_results_all.csv (JÁ EXISTENTE)              │
│     ├── Filtra compostos nomeados → RESULTADO VAZIO                  │
│     └── Gráfico de barras top genes × compostos                     │
│                                                                      │
│  7. REDE PPI (STRING PRONTO)                                         │
│     ├── Carrega DEG_AMPK_pathway.csv (JÁ EXISTENTE)                 │
│     ├── Carrega STRING protein.info (JÁ BAIXADO)                    │
│     ├── Carrega STRING protein.links (JÁ BAIXADO)                   │
│     ├── Mapeia 119 genes AMPK → rede                                │
│     └── Gera grafo igraph + figura ggplot2                           │
│                                                                      │
│  8. RELATÓRIO DE AUDITORIA                                           │
│     ├── audit_events_nqo1_pubchem_ppi.csv                           │
│     └── AUDIT_NQO1_PUBCHEM_PPI.md                                   │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

### O que NÃO está no script.R (quebra de reprodutibilidade)

```
┌─────────────────────────────────────────────────────────────────────┐
│               AUSÊNCIAS CRÍTICAS NO PIPELINE                          │
├─────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  ❌ Download do GSE116250 via GEOquery                                │
│  ❌ Processamento com limma (lmFit, eBayes, topTable)                │
│  ❌ Cálculo de DEG (tabela completa)                                  │
│  ❌ Filtro de genes da via AMPK a partir do DEG completo              │
│  ❌ GSEA (fgsea)                                                     │
│  ❌ GSVA                                                            │
│  ❌ Consulta ao PubChem via API                                       │
│  ❌ Download dos dados do STRING                                      │
│  ❌ Geração do DGB_results_NQO1.xlsx                                 │
│  ❌ Heatmap da via AMPK                                              │
│  ❌ Boxplot GSVA                                                     │
│  ❌ Volcano plot do DEG                                              │
│                                                                      │
│  Todos estes outputs existem nos diretórios mas foram gerados        │
│  por SCRIPT(S) NÃO COMPARTILHADO(S).                                 │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

---

## ETAPA 2, AUDITORIA LINHA POR LINHA DO script.R

### Bloco 1: Cabeçalho e metadados (linhas 1-59)
- **O que faz**: Documenta estrutura de diretórios e instruções de uso
- **Problema**: A estrutura declarada menciona `output/tables/ampk/`, `output/tables/pubchem/`, `output/tables/nqo1/` mas **não menciona** `output/tables/ppi/`, que o código efetivamente cria e usa
- **Gravidade**: BAIXA (inconsistência de documentação)

### Bloco 2: Configuração de diretórios (linhas 60-81)
- **O que faz**: Detecta raiz do projeto, define estrutura de diretórios
- **Problema**: `tables_ppi` está definido no código (linha 71) mas NÃO está documentado no cabeçalho
- **Gravidade**: BAIXA

### Bloco 3: Sistema de auditoria (linhas 83-96)
- **O que faz**: Cria `AUDIT_EVENTS` como data.frame global, função `audit_log()` com `<<-`
- **Problema**: Uso de `<<-` (superassignment) é correto aqui, mas o data.frame cresce com `rbind()` a cada evento, O(n²) ineficiente, mas com ~15 eventos é irrelevante
- **Gravidade**: BAIXA

### Bloco 4: Carregamento de pacotes (linhas 98-106)
- **O que faz**: `require_pkg()` verifica e carrega dplyr, ggplot2, stringr, openxlsx, igraph
- **OK**: Sem problemas

### Bloco 5: Funções utilitárias (linhas 108-139)
- **`find_first_existing`**: Busca arquivo em múltiplos paths. OK.
- **`save_plot`**: Salva com ggsave. OK.
- **`read_csv_flexible`**: Lê CSV com trimws nos nomes. OK.
- **`standardize_names`**: Converte espaços/hífens/pontos para underscore. OK.
- **⚠️ `cleanup_unneeded_outputs`**: Lista hardcoded de arquivos para deletar. Remove `DGB_results_NQO1.csv` e `DGB_results_NQO1_clean.csv` mas o script está gerando `DGB_results_NQO1_from_xlsx_clean.csv`, nomes diferentes. OK (não conflita).

### Bloco 6: `decimal_to_numeric` (linhas 141-149)
- **O que faz**: Converte números com vírgula decimal (formato europeu) para numeric
- **Lógica**: Se tem vírgula, remove todos os pontos (separadores de milhar) e troca vírgula por ponto
- **⚠️ Problema**: Se um número tem ponto como decimal E vírgula como milhar (ex: "1.234,56"), funciona. Mas se o formato for "1,234.56" (vírgula como milhar, ponto como decimal, formato anglo-americano), a função removeria a vírgula incorretamente. Isso depende do locale do Excel que gerou o XLSX.
- **Gravidade**: BAIXA (funciona para o XLSX atual, mas frágil)

### Bloco 7: `q_to_neglog10` (linhas 151-154)
- **O que faz**: Calcula -log10(q) com cap em 50
- **OK**: Correta. Valores q=0 são capados em 10⁻⁵⁰ → -log10 = 50

### Bloco 8: `extract_dose_um` (linhas 156-167)
- **O que faz**: Extrai valor numérico e unidade (nM, µM, mM) de string de dose, converte para µM
- **Regex**: `\\d+(?:\\.\\d+)?`, captura números com ponto decimal
- **⚠️ Problema**: Se a dose vier com vírgula decimal (ex: "0,1 µM"), o regex `\\d+(?:\\.\\d+)?` NÃO captura porque espera ponto. Mas `str_replace_all(dose, ",", ".")` é aplicado antes.
- **OK**: Funciona para o formato atual, mas a ordem das operações é frágil

### Bloco 9: `read_nqo1_xlsx` (linhas 169-218), FUNÇÃO CENTRAL
- **O que faz**: Lê todas as abas do XLSX, padroniza colunas, classifica direção, calcula score
- **Verificações**:
  - ✅ Exige colunas Drug_Name, p_value, q_value, Fold_Change
  - ✅ Detecta fonte (CMAP/LINCS/CREEDS) pelo nome da aba
  - ✅ Detecta expected_direction (up/down) pelo sufixo da aba
  - ⚠️ `expected_direction` é calculado mas **NUNCA USADO** depois, código morto
  - ⚠️ Se `Specificity` não existe, define como NA, mas Specificity **nunca é usado** depois
  - ✅ `regulation` é definido por `Fold_Change > 0` → "Upregula", `< 0` → "Downregula"
  - ✅ `score = neg_log10_q_capped * abs(Fold_Change)`
  - ⚠️ **`score` usa `neg_log10_q_capped` (que é baseado em `q_value_for_plot = pmax(q_value, 1e-50)`)**, mas `q_to_neglog10` usa `pmax(q_value, 10^-cap)` onde cap=50. São consistentes.

### Bloco 10: `rank_nqo1_drugs` (linhas 220-243), FUNÇÃO CRÍTICA
- **O que faz**: Agrupa por Drug_Name, calcula estatísticas resumo
- **⚠️ PROBLEMA**: Agrupamento é **CASE SENSITIVE**. "Doxorubicin" ≠ "doxorubicin"
- **Evidência nos outputs**:
  ```
  Doxorubicin (CREEDS):     n_tests=10, direction="Downregula NQO1"
  doxorubicin (LINCS L1000): n_tests=5,  direction="Downregula NQO1"
  ```
- **direction** usa `median_fc`, não `mean_fc`. Correto pela lógica, mas manumycin-a tem mean_fc negativo (-0,0097) e median_fc positivo (+0,746) → classificado como "Upregula", contra-intuitivo.
- **Gravidade**: MÉDIA

### Bloco 11: `plot_volcano` (linhas 245-275)
- **O que faz**: Volcano plot com pontos coloridos por regulação
- **q_was_zero**: trichostatin A tem q=0 → formato diferente (triângulo)
- **Top labels**: 12 fármacos com menor q-value e maior |Fold_Change|
- **OK**: Código correto

### Bloco 12: `plot_dose_response` (linhas 277-344)
- **O que faz**: Dose-response facetado por fármaco
- **Lógica complexa de seleção**:
  1. Filtra drogas com ≥2 doses numéricas → ordena por q-value e score → top 9
  2. Força inclusão de Doxorubicin (detectada por string "doxo")
  3. Fallback: se nenhuma droga tem ≥2 doses, usa top por q-value
- **OK**: Código robusto, com fallbacks

### Bloco 13: `load_pubchem` (linhas 346-356)
- **O que faz**: Carrega CSV PubChem pré-existente
- **OK**: Se não encontra, retorna data.frame vazio

### Bloco 14: `plot_pubchem_counts` (linhas 358-375)
- **O que faz**: Gráfico de barras top 20 genes × nº de compostos
- **OK**: Se pubchem_df vazio, retorna NULL

### Bloco 15: `plot_ppi_network` (linhas 377-462)
- **O que faz**: Constrói rede PPI a partir de DEG + STRING
- **Fluxo**:
  1. Lê DEG_AMPK_pathway.csv
  2. Mapeia genes → STRING protein IDs via protein.info
  3. Filtra links com combined_score (do STRING)
  4. Constrói grafo igraph
  5. Plota com ggplot2 (geom_segment + geom_point + geom_text)
- **⚠️** O grafo usa **todos os links**, não apenas os de alta confiança (>700). O `combined_score` é usado como espessura da linha (linewidth), mas não como filtro. O filtro min700 foi aplicado no download do arquivo (`.min700.txt.gz`), não no código.
- **OK**: Correto, o filtro já está no nome do arquivo

### Bloco 16: `main()` (linhas 464-530)
- **Sequência**: cleanup → ler XLSX → ranking → salvar CSVs → volcano → dose-response → PubChem → PPI → relatório
- **`doxo` filter (linha 491)**: Busca "doxo|doxorubic|doxorubicin" no ranking
  - ⚠️ Encontra apenas "Doxorubicin" (case sensitive match no ranking). "doxorubicin" (minúsculo) NÃO é capturado por este filtro porque `str_to_lower` é aplicado a `Drug_Name` mas o texto buscado é "doxo" que casa com ambos.
  - ✅ Na verdade, lendo com mais cuidado: `str_detect(str_to_lower(Drug_Name), "doxo|doxorubic|doxorubicin")`, isso captura ambas as grafias. Mas o ranking tem as duas entradas separadas, então o `doxo` CSV terá 2 linhas.

---

## ETAPA 3, AUDITORIA METODOLÓGICA (REPRODUTIBILIDADE)

### Matriz de reprodutibilidade

| Método declarado no paper | Presente no script.R? | Output existe? | REPRODUTÍVEL? |
|---------------------------|----------------------|----------------|---------------|
| GEOquery (download GSE116250) | ❌ NÃO | GSE116250_* em data/raw/ | ❌ Script externo |
| limma (DEG) | ❌ NÃO | DEG_full_table.csv, DEG_AMPK_pathway.csv | ❌ Script externo |
| GSEA (fgsea) | ❌ NÃO | GSEA_AMPK_results.csv, GSEA_KEGG_all_pathways.csv | ❌ Script externo |
| GSVA | ❌ NÃO | GSVA_AMPK_scores_per_sample.csv | ❌ Script externo |
| STRING (rede PPI) | ✅ Parcial | protein.links, protein.info em data/raw/ | ⚠️ Lê dados prontos |
| Drug Gene Budger (DGB) | ❌ NÃO | DGB_results_NQO1.xlsx | ❌ Ferramenta externa |
| PubChem | ❌ NÃO | pubchem_results_all.csv | ❌ Script externo |
| Processamento NQO1 | ✅ SIM | DGB_results_NQO1_from_xlsx_clean.csv | ✅ |
| Ranking fármacos | ✅ SIM | ranking_drugs_NQO1.csv | ✅ |
| Volcano plot | ✅ SIM | volcano_NQO1_drugs.png | ✅ |
| Dose-response | ✅ SIM | dose_response_NQO1_drugs.png | ✅ |
| Rede PPI (igraph/ggplot2) | ✅ SIM | PPI_AMPK_network.png | ✅ |

**Conclusão**: O script.R cobre apenas a **porção final do pipeline** (pós-processamento e visualização). As etapas fundamentais de bioinformática (DEG, GSEA, GSVA, consulta PubChem, geração do XLSX do DGB) foram executadas em **scripts não compartilhados**. A metodologia descrita no paper **não pode ser reproduzida** apenas com os arquivos disponibilizados.

---

## ETAPA 4, AUDITORIA DOS RESULTADOS (CRUZAMENTO PAPER × OUTPUTS)

### Tabela de todas as afirmações numéricas do paper × dados reais

| Afirmação no Paper | Valor no Paper | Valor Real nos Outputs | Status |
|-------------------|----------------|------------------------|--------|
| q-value da doxorrubicina | 6,32×10⁻⁶ | 2,21×10⁻²⁸ (CREEDS) / 2,50×10⁻⁹ (LINCS) | 🔴 FALSO |
| Doxorrubicina "indutor direto" | Upregulation | Downregulation (mean_fc = −0,315 / −2,565) | 🔴 INVERSÃO DE SINAL |
| Geldanamicina "elevação transcricional" | Upregulation | Downregulation (mean_fc = −0,275) | 🔴 INVERSÃO DE SINAL |
| Wortmannin "elevação transcricional" | Upregulation | Downregulation (mean_fc = −0,159) | 🔴 INVERSÃO DE SINAL |
| Monensina "potente capacidade regulatória" | Mencionada | NÃO CONSTA NOS DADOS | 🔴 INEXISTENTE |
| NQO1 como alvo da "triagem da via AMPK" | Implícito | NQO1 NÃO está na via AMPK (KEGG) | 🔴 PREMISSA FALSA |
| GSEA da via AMPK | Não mencionado | p=0,125, padj=0,419 (NÃO SIGNIFICATIVO) | 🟡 OMISSÃO |
| "3523 linhas válidas" no DGB | Não mencionado | 3523 (confirmado) | ✅ |
| "1281 fármacos únicos" | Não mencionado | 1281 (confirmado) | ✅ |
| Top 1 do ranking = trichostatin A | Não mencionado | trichostatin A (confirmado) | ✅ |
| Alvocidib "potente repressor" | Downregulation | Downregulation (mean_fc = −1,485) | ✅ |
| Doxorrubicina a 10 µM | 10 µM | Confirmado nos dados (LINCS, 10 µM) | ✅ |

---

## ETAPA 5, AUDITORIA BIOLÓGICA

### 5.1 NQO1 pertence à via AMPK?

**NÃO.** Evidências:

1. `AMPK_pathway_genes_KEGG.csv` contém 119 genes da via AMPK (hsa04152). **NQO1 NÃO está na lista.**
2. `DEG_AMPK_pathway.csv` contém apenas genes da via AMPK. **NQO1 NÃO aparece.**
3. NQO1 só aparece em `DEG_full_table.csv` (tabela completa de todos os genes do microarray), onde é o 2º gene mais significativo (padj = 6,86×10⁻¹³), mas **isso não tem relação com a via AMPK**.

### 5.2 A narrativa "triagem AMPK → NQO1" se sustenta?

**NÃO.** A narrativa do paper constrói o seguinte encadeamento:
> "triagem da via AMPK" → "genes diferencialmente expressos" → "dentre os alvos, destaca-se NQO1"

Isso é uma **construção narrativa artificial**. O NQO1 foi selecionado da tabela DEG completa (todos os ~42.000 genes), não do subconjunto AMPK. O paper cria uma falsa conexão metodológica entre a via AMPK e o NQO1.

### 5.3 Evidência para conectar AMPK com NQO1?

**Muito tênue.** Existe literatura mostrando que AMPK pode regular NRF2, que por sua vez regula NQO1. Mas:
- Isso não é mediado diretamente pela via canônica AMPK do KEGG
- O paper não cita nem desenvolve essa conexão mecanicista
- A justificativa se limita a "NQO1 guarda estrita relação com mecanismos de cardiotoxicidade" (referência [7] sobre doxorrubicina, não sobre AMPK)

### 5.4 Extrapolação biológica presente?

**SIM, GRAVE.** O paper extrapola de:
- Dados de linhagens tumorais (MCF7, PC3, VCAP) para cardiomiócitos humanos
- Associação estatística in vitro para "mecanismo compensatório homeostático"
- Dados de perturbação aguda (6h-24h) para "farmacovigilância" e "cuidados assistenciais de precisão"

### 5.5 Cherry-picking?

**SIM.** O paper:
- Seleciona NQO1 dentre ~42.000 genes sem justificar claramente que NÃO veio da via AMPK
- Omite que o GSEA da via AMPK não foi significativo
- Cita apenas fármacos que "confirmam" a narrativa (doxorrubicina como indutor), ignorando que ela é DOWNREGULADORA

### 5.6 Overinterpretation?

**SIM, GENERALIZADA.** Exemplos:
- "fornece bases teóricas sólidas" → dados in vitro de linhagens tumorais
- "subsidiam a farmacovigilância" → sem validação clínica
- "enfermagem de precisão no monitoramento do risco cardiometabólico" → extrapolação translacional extrema

---

## ETAPA 6, AUDITORIA ESTATÍSTICA

### 6.1 Múltiplas comparações
- **GSEA**: 163 vias testadas. Correção FDR aplicada (coluna padj). ✅ OK.
- **DEG**: ~42.000 genes. Correção FDR aplicada (coluna padj). ✅ OK.
- **DGB**: 3523 testes fármaco-condição. q-values vêm da plataforma DGB (não está claro qual correção). ⚠️

### 6.2 Interpretação de q-value
- O paper cita "q-value = 6,32×10⁻⁶" para doxorrubicina. Este valor **não existe**. O q-value real é 2,21×10⁻²⁸. Isso não é um erro de interpretação, é um **dado incorreto**.

### 6.3 Interpretação de Fold Change
- O paper interpreta Fold Change positivo como "indução". Mas o Fold Change no DGB representa a mudança na expressão de NQO1 após tratamento. Para doxorrubicina, é **negativo** (−0,315 e −2,565), indicando **repressão**, não indução.

### 6.4 Volcano plot
- O volcano plot gerado pelo código está correto. Mostra claramente doxorrubicina no lado esquerdo (Fold Change negativo). A legenda diz "Downregula NQO1". **O paper contradiz diretamente a figura.**

### 6.5 Score de ranking
- `score = -log10(q) × |Fold_Change|`. É uma heurística razoável, mas não tem fundamento estatístico formal. Deveria ser justificado.

### 6.6 Desbalanceamento amostral GSVA
- CTRL: n=14, DCM: n=37. Desbalanceamento de 2,6×. Testes estatísticos paramétricos podem ser afetados.

---

## ETAPA 7, AUDITORIA FARMACOGENÔMICA

### 7.1 Direção da regulação

| Fármaco | Paper diz | Dados mostram | Correto? |
|---------|-----------|---------------|----------|
| Doxorrubicina | ⬆️ Upregulation | ⬇️ Downregulation | 🔴 INVERSÃO |
| Geldanamicina | ⬆️ Upregulation | ⬇️ Downregulation | 🔴 INVERSÃO |
| Wortmannin | ⬆️ Upregulation | ⬇️ Downregulation | 🔴 INVERSÃO |
| Alvocidib | ⬇️ Downregulation | ⬇️ Downregulation | ✅ |
| Monensina | Mencionada | INEXISTENTE | 🔴 FANTASMA |

### 7.2 Compostos inexistentes
- **Monensina**: Zero registros no ranking. Busca textual: nenhum resultado.

### 7.3 Ranking do paper vs código
O paper menciona "cinco primeiros colocados no ranking molecular". O ranking real:
1. trichostatin A (CMAP), q=0, Up
2. Doxorubicin (CREEDS), q=2,21e-28, Down
3. manumycin-a (LINCS), q=3,62e-20, Up
4. parthenolide (LINCS), q=4,42e-19, Up
5. panobinostat (LINCS), q=6,60e-18, Down

O paper **não menciona trichostatin A, manumycin-a, parthenolide nem panobinostat**, foca apenas em doxorrubicina, geldanamicina, wortmannin e monensina, que **não estão entre os top 5** (geldanamicina está na posição ~48, wortmannin ~60, monensina não existe).

---

## ETAPA 8, AUDITORIA DAS FIGURAS

### 8.1 `volcano_NQO1_drugs.png`
- **Corresponde aos dados?** ✅ Sim
- **Corresponde ao código?** ✅ Sim
- **Eixos corretos?** ✅ FC × -log10(q)
- **Cores:** ✅ Verde = up, Laranja = down
- **Interpretação no paper:** 🔴 O paper descreve o volcano citando compostos que aparecem no lado oposto ao descrito

### 8.2 `dose_response_NQO1_drugs.png`
- **Corresponde aos dados?** ✅ Sim
- **Subtítulo:** "Doxorubicin incluída a partir dos registros LINCS com 10 µM", ✅
- **OK**

### 8.3 `pubchem_top_genes_compound_counts.png`
- **Dados de origem:** `pubchem_results_all.csv` (1887 linhas, 1874 compostos sem nome)
- **Gráfico gerado:** Top genes com mais compostos
- **OK**, mas `pubchem_results_with_names.csv` está vazio

### 8.4 `PPI_AMPK_network.png`
- **119 nós, 1874 arestas** ✅
- **Cor:** log2FC (verde = up, laranja = down) ✅
- **OK**

### 8.5 Figuras AMPK (não geradas pelo script.R)
- `volcano_plot.png`, `heatmap_AMPK_pathway.png`, `heatmap_top_DEGs.png`, `GSEA_AMPK_enrichment_plot.png`, `GSVA_AMPK_boxplot.png`
- Geradas por script externo. Não auditáveis completamente.

---

## ETAPA 9, AUDITORIA TEXTUAL (FRASE POR FRASE)

### Resumo

| # | Frase | Classificação | Motivo |
|---|-------|--------------|--------|
| 1 | "triagem transcriptômica in silico da via AMPK utilizando dados do dataset GSE116250" | PARCIALMENTE CORRETA | O GSE116250 foi usado, mas a triagem não foi "da via AMPK", NQO1 veio do DEG completo |
| 2 | "gene NQO1 foi selecionado devido à sua relevância biológica na citoproteção e homeostase redox" | CORRETA | NQO1 tem esse papel documentado |
| 3 | "doxorrubicina atua como um importante modulador transcriptômico positivo" | **FALSA** | Dados mostram DOWNregulation (mean_fc = −0,315 e −2,565) |
| 4 | "q-value = 6,32e-6" | **FALSA** | Valor real: 2,21×10⁻²⁸. Nenhum registro tem q=6,32e-6 |
| 5 | "sugerindo o acionamento de uma resposta adaptativa ou compensatória" | **ESPECULATIVA** | Inferência mecanística sem experimento funcional; direção oposta aos dados |
| 6 | "Geldanamicina e Monensina" como potentes reguladores | **FALSA** (Geldanamicina: direção errada) / **FALSA** (Monensina: inexistente) | Geldanamicina é downreguladora; Monensina não consta |
| 7 | "fornece bases teóricas sólidas que subsidiam a farmacovigilância" | **EXAGERADA** | Dados in vitro de linhagens tumorais não subsidiam farmacovigilância sem validação |

### Introdução

| # | Frase | Classificação |
|---|-------|--------------|
| 8 | "CMD consiste em uma das principais etiologias" | CORRETA |
| 9 | "AMPK funciona como um sensor energético celular mestre" | CORRETA |
| 10 | "triagem sistemática baseada em biologia de sistemas da via AMPK" | **NÃO SUPORTADA**, a triagem que encontrou NQO1 não foi restrita à via AMPK |
| 11 | "dentre os alvos desregulados identificados no rastreio da via, destaca-se o gene NQO1" | **FALSA**, NQO1 NÃO foi identificado no rastreio da via AMPK |

### Métodos

| # | Frase | Classificação |
|---|-------|--------------|
| 12 | "análise estatística de expressão diferencial [...] pacotes GEOquery e limma" | **IRREPRODUZÍVEL**, Não está no script.R |
| 13 | "ajuste de p-valor pela taxa de falsa descoberta (FDR, significância para padj < 0,05)" | PARCIALMENTE CORRETA, existe nos outputs, mas o código não está no script.R |
| 14 | "filtrados no STRING utilizando um escore de confiança restritivo mínimo de 700" | CORRETA, o arquivo é `.min700.txt.gz` |
| 15 | "modelagem topológica da rede proteína-proteína (PPI) foi executada via pacote igraph" | CORRETA |
| 16 | "gene NQO1 isolado na triagem foi submetido à varredura de perturbações transcricionais no DGB" | PARCIALMENTE CORRETA, o DGB foi usado, mas a "triagem" não isolou NQO1 da via AMPK |
| 17 | "computou-se um escore analítico customizado" | CORRETA, `score = -log10(q) × |FC|` |

### Resultados e Discussão

| # | Frase | Classificação |
|---|-------|--------------|
| 18 | "triagem transcriptômica inicial sinalizou alterações marcantes na expressão de fatores moleculares constituintes da via da AMPK" | PARCIALMENTE CORRETA, há DEGs na via, mas GSEA não significativo |
| 19 | "gene NQO1 despontou como um componente periférico estratégico" | **INCONSISTENTE**, NQO1 não está na via AMPK, não pode "despontar" dela |
| 20 | "sua desregulação em pacientes com CMD correlaciona-se com o estresse miocárdico" | **ESPECULATIVA**, não há análise de correlação NQO1 × estresse |
| 21 | "capacidade biológica intrínseca de ligar o controle metabólico da AMPK ao tamponamento de radicais livres" | **ESPECULATIVA**, inferência sem dados |
| 22 | "Geldanamicina e o Wortmannin induzem expressiva elevação transcricional sobre o NQO1" | **FALSA**, ambos são DOWNREGULADORES |
| 23 | "doxorrubicina demonstrou associação estatisticamente significante como indutor direto da transcrição de NQO1 (q-value = 6,32 10 e-6)" | **FALSA DUPLA**, direção errada E q-value falso |
| 24 | "registrando maior estabilidade regulatória na concentração de 10" | **INCOMPLETA**, 10 µM, faltou unidade |
| 25 | "indução transcriptômica sugere a ativação de uma via adaptativa miocárdica" | **FALSA**, a direção é de repressão, não indução |

### Conclusões

| # | Frase | Classificação |
|---|-------|--------------|
| 26 | "triagem transcriptômica em ambiente computacional permitiu selecionar o gene NQO1" | PARCIALMENTE CORRETA |
| 27 | "NQO1 como um elemento metabólico periférico e protetor relevante na fisiopatologia da CMD" | **ESPECULATIVA** |
| 28 | "doxorrubicina atua como um regulador transcricional positivo da expressão de NQO1" | **FALSA**, é downreguladora |
| 29 | "indicando o acionamento de um mecanismo compensatório homeostático celular de proteção redox" | **FALSA**, premissa errada (down, não up) |
| 30 | "fornecem insumos conceituais para balizar ações preventivas de farmacovigilância" | **EXAGERADA** |

---

## ETAPA 10, FALÁCIAS E VIESES IDENTIFICADOS

### 10.1 Falácia ecológica
- **Dados de linhagens tumorais** (MCF7 = câncer de mama, PC3 = próstata, VCAP = próstata) são usados para inferir comportamento de **cardiomiócitos humanos**. O microambiente celular, perfil epigenético e vias de sinalização são completamente diferentes.

### 10.2 Falácia causal
- Associação estatística fármaco → expressão de NQO1 é tratada como **relação causal direta** ("atua como modulador", "indução direta", "acionamento de mecanismo compensatório").

### 10.3 Viés de confirmação
- O paper seleciona apenas fármacos que "confirmam" a narrativa desejada (doxorrubicina, geldanamicina, wortmannin, monensina), ignorando os top 5 reais do ranking (trichostatin A, manumycin-a, parthenolide, panobinostat).

### 10.4 Cherry-picking
- NQO1 selecionado do DEG completo, mas apresentado como fruto da "triagem da via AMPK"
- GSEA não significativo da via AMPK completamente omitido
- Apenas fármacos familiares mencionados; top 1 (trichostatin A) ignorado

### 10.5 Spin científico
- Dados que mostram DOWNregulation são apresentados como UPregulation em TODAS as instâncias principais
- Resultado negativo (GSEA AMPK p=0,125) é omitido
- q-value falso (6,32e-6) é apresentado como evidência de significância

### 10.6 HARKing (Hypothesizing After Results are Known)
- A conexão AMPK → NQO1 parece ter sido construída post hoc para dar uma narrativa metodológica ao achado do NQO1 no DEG completo

### 10.7 Exagero translacional
- Dados puramente in silico → "farmacovigilância", "enfermagem de precisão", "cuidados assistenciais", "cardiologia translacional moderna"

---

## ETAPA 11, AUDITORIA DAS REFERÊNCIAS

| Ref | Tema | Suporta? | Observação |
|-----|------|----------|------------|
| [1] BOZKURT, 2021, CMD review | ✅ Parcialmente | Contexto geral de CMD |
| [2] RAMACCINI, 2021, Mitochondrial function in DCM | ✅ Parcialmente | Mitocôndria em CMD |
| [3] MARTINS, 2022, Oxidative stress in cardiac remodeling | ✅ Parcialmente | Estresse oxidativo |
| [4] HERZIG & SHAW, 2018, AMPK review | ✅ | AMPK, mas não conecta com NQO1 |
| [5] ROSS & SIEGEL, 2017, NQO1 functions | ✅ | Papel do NQO1 |
| [6] HODOS, 2016, In silico drug repurposing | ✅ | Métodos in silico |
| [7] WALLACE, 2003, Doxorubicin cardiotoxicity | ✅ Parcialmente | Doxorrubicina, mas artigo de 2003 |
| [8] WANG, 2021, Drug Gene Budger | ✅ | Ferramenta DGB |
| [9] EKINS, 2019, ML for drug discovery | ⚠️ Fraca | Não diretamente relacionado |
| [10] ZANGIROLAMI-RAIMUNDO, 2018, Estudos transversais | ⚠️ Fraca | Metodologia geral, não específica |
| [11] COSTA, 2026, Resultados de bioinformática IFF | ⚠️ Não verificável | Publicação interna |
| [12] LOPES-JÚNIOR, 2021, Personalized nursing | ⚠️ Tangencial | Enfermagem de precisão |
| [13] SOUZA & PARDINI, 2016, Carbono pirolítico | 🔴 **TOTALMENTE DESCONEXA** | **Engenharia de materiais. Zero relação com o tema.** |

**Referência [13] é uma contaminação de outro trabalho.** "Identificação do Tipo de Carbono Pirolítico por Meio da Técnica do Ângulo de Extinção" é sobre engenharia de materiais no ITA/IAE. Foi colada indevidamente.

---

## ETAPA 12, CONSISTÊNCIA INTERNA

### Inconsistências entre seções do paper

| Elemento A | Elemento B | Inconsistência |
|-----------|-----------|----------------|
| Resumo: "q-value = 6,32e-6" | Dados: q-value = 2,21e-28 | 🔴 Valor diferente |
| Resumo: "modulador positivo" | Dados: Downregulation | 🔴 Sinal invertido |
| Introdução: "triagem da via AMPK" | Dados: NQO1 não está na via | 🔴 Premissa falsa |
| Métodos: "GEOquery e limma" | script.R: não contém | 🔴 Irreproduzível |
| Resultados: "elevação com Geldanamicina e Wortmannin" | Dados: Downregulation | 🔴 Sinal invertido |
| Conclusão: "regulador transcricional positivo" | Dados: Downregulation | 🔴 Sinal invertido |
| Paper: NQO1 "componente periférico" da via | Dados: NQO1 FORA da via | 🔴 Inconsistência |

---

## ETAPA 13, CORREÇÕES SUGERIDAS

### Correção 1: q-value da doxorrubicina

**Problema**: Valor 6,32×10⁻⁶ não existe
**Impacto**: CRÍTICO
**Correção**: Usar o valor real

> **Texto atual**: "q-value = 6,32×10⁻⁶"
> **Texto corrigido**: "q-value = 2,21×10⁻²⁸ (CREEDS) e 2,50×10⁻⁹ (LINCS L1000)"

### Correção 2: Direção da regulação da doxorrubicina

**Problema**: Paper diz "upregulation", dados mostram "downregulation"
**Impacto**: CRÍTICO, Invalida toda a discussão e conclusão
**Novo texto**:
> "A doxorrubicina demonstrou associação estatisticamente significante como **repressor transcricional** de NQO1 (Fold Change médio = −0,315 no CREEDS e −2,565 no LINCS L1000; q-value mínimo = 2,21×10⁻²⁸), indicando que esta antraciclina **suprime** a expressão de NQO1 em modelos celulares, o que pode exacerbar a vulnerabilidade ao estresse oxidativo."

### Correção 3: Geldanamicina e Wortmannin

**Problema**: Descritos como "elevação transcricional"
**Impacto**: CRÍTICO
**Novo texto**:
> "Compostos como a Geldanamicina (Fold Change médio = −0,275; direção: downregulação de NQO1) e Wortmannin (Fold Change médio = −0,159; direção: downregulação de NQO1) operam como **repressores enzimáticos** de NQO1."

### Correção 4: Monensina

**Problema**: Citada mas inexistente nos dados
**Impacto**: CRÍTICO
**Correção**: **Remover completamente a menção à Monensina.**

### Correção 5: Conexão AMPK → NQO1

**Problema**: NQO1 não está na via AMPK
**Impacto**: ALTO
**Novo texto (Introdução)**:
> "A análise de expressão diferencial global (DEG) do transcriptoma completo identificou alterações significativas em múltiplos genes. Paralelamente, genes da via AMPK também apresentaram desregulação, embora a via como um todo não tenha mostrado enriquecimento estatisticamente significativo no GSEA (p = 0,125, padj = 0,419). Dentre os genes mais diferencialmente expressos no transcriptoma completo, destacou-se o NQO1 (padj = 6,86×10⁻¹³, log2FC = −1,686), selecionado para análise farmacogenômica devido à sua relevância na citoproteção redox."

### Correção 6: Referência [13]

**Problema**: Totalmente desconexa
**Impacto**: ALTO
**Correção**: **Remover completamente.** Substituir por referência pertinente ou remover a citação.

### Correção 7: Conclusão

**Problema**: Baseada em premissas falsas
**Impacto**: CRÍTICO
**Novo texto**:
> "A análise transcriptômica revelou que o gene NQO1 está significativamente **reprimido** (log2FC = −1,686, padj = 6,86×10⁻¹³) no miocárdio de pacientes com CMD. A modelagem farmacogenômica indicou que a doxorrubicina atua como **repressor transcricional** de NQO1 em linhagens celulares, sugerindo que este quimioterápico pode agravar o déficit de NQO1 já presente na CMD. Estes achados in silico, embora requeiram validação experimental em cardiomiócitos primários, apontam para um potencial mecanismo de exacerbação da cardiotoxicidade por antraciclinas mediado pela supressão de defesas antioxidantes. A transposição destes achados para a prática clínica requer estudos funcionais e validação em coortes independentes."

---

## ETAPA 14, AVALIAÇÃO DE CREDIBILIDADE CIENTÍFICA

| Critério | Nota (0-10) | Justificativa |
|----------|-------------|---------------|
| Originalidade | 5 | Análise farmacogenômica de NQO1 é relevante, mas abordagem é padrão |
| Metodologia | 2 | Script principal não contém etapas-chave; irreproduzível |
| Bioinformática | 3 | Pipeline fragmentado; outputs de scripts externos sem código-fonte |
| Reprodutibilidade | 1 | Script.R cobre apenas ~20% do pipeline descrito |
| Estatística | 3 | q-value falso no paper; correções existem nos outputs mas não no texto |
| Interpretação biológica | 1 | Direção da regulação invertida para 3 fármacos; premissa AMPK→NQO1 falsa |
| Discussão | 1 | Baseada em dados incorretos; extrapolação translacional excessiva |
| Conclusões | 1 | Invalidada por inversão de sinal e q-value falso |
| Farmacogenômica | 2 | Ranking ignorado; fármacos citados não conferem |
| Qualidade das figuras | 6 | Figuras geradas estão corretas; contraditas pelo texto |
| Consistência interna | 1 | Múltiplas contradições resumo↔dados, métodos↔código |
| Potencial de publicação | 1 | Não publicável no estado atual |

**Média geral: 2,3/10**

---

## ETAPA 15, PARECER FINAL (FORMATO NATURE REVIEWER)

### RECOMMENDATION: ❌ REJECT, MAJOR REVISION REQUIRED BEFORE RESUBMISSION

---

### PRINCIPAIS FORÇAS
1. Pergunta de pesquisa relevante: cardiotoxicidade por antraciclinas é um problema clínico real
2. Uso de bases públicas (GEO, STRING) com potencial de reprodutibilidade
3. Integração de múltiplas plataformas farmacogenômicas (CMAP, LINCS, CREEDS)
4. Código da parte final do pipeline está organizado e documentado

---

### ERROS FATAIS (IMPEDEM PUBLICAÇÃO)

**F1. q-value falso no manuscrito**
O valor 6,32×10⁻⁶ atribuído à doxorrubicina não corresponde a nenhum registro nos dados. O valor real é 2,21×10⁻²⁸. Isso configura **dado incorreto no manuscrito**.

**F2. Direção da regulação invertida para 3 fármacos**
Doxorrubicina, Geldanamicina e Wortmannin são DOWNREGULADORES de NQO1 nos dados, mas o manuscrito afirma o contrário. A tese central do paper (resposta adaptativa/compensatória mediada por upregulation de NQO1) é **invalidada pelos próprios dados**.

**F3. Monensina inexistente**
O composto é citado como achado mas não consta em nenhum dos 1281 fármacos do ranking. Isso sugere **inclusão de dado não verificado**.

**F4. Premissa AMPK→NQO1 artificial**
NQO1 não pertence à via AMPK do KEGG. O paper constrói uma falsa conexão metodológica. A análise que identificou NQO1 foi do transcriptoma completo, não da via AMPK.

**F5. Irreprodutibilidade do pipeline**
O script.R cobre apenas pós-processamento. As etapas de GEOquery, limma, GSEA, GSVA e consulta PubChem estão ausentes. O estudo não é reproduzível com os arquivos fornecidos.

---

### ERROS GRAVES (COMPROMETEM A CREDIBILIDADE)

**G1.** GSEA da via AMPK não significativo (p=0,125, padj=0,419), omitido do paper

**G2.** Referência [13] é de engenharia de materiais (carbono pirolítico), completamente desconexa

**G3.** Afirmações de "farmacovigilância", "enfermagem de precisão" e "cardiologia translacional" são exageros translacionais sem suporte

**G4.** Ranking top 5 real (trichostatin A, doxorubicin, manumycin-a, parthenolide, panobinostat) não é discutido

**G5.** `pubchem_results_with_names.csv` vazio, mapeamento PubChem não gerou resultados nomeados

---

### O QUE PRECISA SER CORRIGIDO OBRIGATORIAMENTE

1. ✅ Substituir q-value 6,32e-6 pelo valor real (2,21e-28)
2. ✅ Corrigir direção da regulação: doxorrubicina é DOWNREGULADORA
3. ✅ Corrigir direção de Geldanamicina e Wortmannin: DOWNREGULADORES
4. ✅ Remover Monensina
5. ✅ Corrigir narrativa AMPK→NQO1: NQO1 veio do DEG completo
6. ✅ Reportar GSEA AMPK não significativo
7. ✅ Remover referência [13]
8. ✅ Incluir script completo da análise DEG
9. ✅ Reescrever Discussão e Conclusão com base nos dados reais
10. ✅ Suprimir ou qualificar fortemente alegações translacionais

---

### O QUE PODE SER REMOVIDO
- Menção à Monensina
- Referência [13]
- Afirmações sobre "enfermagem de precisão" e "farmacovigilância" (ou movê-las para "perspectivas futuras")
- Narrativa de "triagem da via AMPK" como origem do NQO1

---

### O QUE DEVE SER REESCRITO
- Resumo (inteiro)
- Resultados e Discussão (inteiro)
- Conclusão (inteira)
- Trecho da Introdução sobre AMPK→NQO1

---

### O QUE DEVE SER REANALISADO
- Unificação de "Doxorubicin" e "doxorubicin" (case sensitive) no ranking
- Justificativa para seleção de fármacos discutidos (por que ignorar trichostatin A, parthenolide, etc?)

---

### O QUE PRECISA DE EXPERIMENTAÇÃO ADICIONAL
- Validação in vitro em cardiomiócitos primários
- Validação em modelos animais de cardiotoxicidade por doxorrubicina
- Análise de polimorfismos de NQO1 em coortes clínicas

---

### QUAIS CONCLUSÕES PERMANECEM VÁLIDAS APÓS A AUDITORIA

1. ✅ NQO1 está significativamente desregulado em CMD (log2FC = −1,686, padj = 6,86e-13)
2. ✅ Doxorrubicina modula a expressão de NQO1 (com direção **oposta** à afirmada)
3. ✅ A rede PPI da via AMPK mostra alterações topológicas em CMD
4. ✅ Existem fármacos com efeitos significativos sobre NQO1 em múltiplas plataformas

---

### QUAIS CONCLUSÕES DEVEM SER REMOVIDAS COMPLETAMENTE

1. ❌ "Doxorrubicina atua como modulador transcriptômico positivo/indutor direto"
2. ❌ "Acionamento de resposta adaptativa/compensatória"
3. ❌ "Geldanamicina e Wortmannin induzem elevação transcricional"
4. ❌ "Monensina como potente regulador"
5. ❌ "Bases teóricas sólidas para farmacovigilância" (no estado atual)
6. ❌ Toda inferência de upregulation como mecanismo protetor

---

## 📋 SUMÁRIO EXECUTIVO

| Categoria | Contagem |
|-----------|----------|
| Erros **CRÍTICOS** (invalidam conclusões) | **4** |
| Erros **ALTOS** (comprometem credibilidade) | **5** |
| Erros **MÉDIOS** (enfraquecem rigor) | **3** |
| Erros **BAIXOS** (cosméticos/técnicos) | **3** |
| Referências inválidas/desconexas | **1** |
| Falácias identificadas | **7** |
| Afirmações FALSAS no texto | **8** |
| Afirmações ESPECULATIVAS/EXAGERADAS | **6** |

---

## 🔴 VEREDITO FINAL

**O manuscrito, na forma atual, NÃO atende aos padrões mínimos de credibilidade científica exigidos para publicação.**

As incongruências entre o texto e os dados são numerosas e graves: valor de q-value inventado, direção da regulação invertida para 3 fármacos, menção a composto inexistente, premissa metodológica artificial e pipeline irreproduzível.

**Recomendação**: Rejeitar na forma atual. Convidar os autores a submeterem versão corrigida após:
1. Correção de TODOS os valores numéricos
2. Reversão da direção da regulação em todo o texto
3. Inclusão do código-fonte completo e reproduzível
4. Reescrita da Discussão e Conclusão com base nos dados reais
5. Remoção de alegações translacionais não suportadas

---

*Relatório gerado em 2026-06-26 por auditoria forense completa de 15 etapas.*
*Todos os achados são baseados em verificação direta dos arquivos fornecidos.*
