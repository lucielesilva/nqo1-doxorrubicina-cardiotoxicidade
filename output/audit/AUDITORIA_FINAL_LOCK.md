# 🔒 AUDITORIA CIENTÍFICA FORENSE DEFINITIVA, LOCK FINAL

**Manuscrito auditado:** `resumo_expandido_luciele_cardiologia_2026_CORRIGIDO.docx`
**Data:** 2026-06-26
**Versão do script:** script.R (corrigido, NÃO re-executado, outputs refletem versão anterior)

---

## ETAPA 1, RECONSTRUÇÃO DO PIPELINE

```
DGB_results_NQO1.xlsx (6 abas: cmap/l1000/creeds × up/down)
    │
    ▼
read_nqo1_xlsx() → 3523 registros
    │  • standardize_names()
    │  • decimal_to_numeric()
    │  • source_database detection
    │  • q_value_for_plot, neg_log10_q_capped, score
    │  • regulation = f(Fold_Change)
    │  • Drug_Name_Norm = str_to_lower()
    ▼
rank_nqo1_drugs() → 1266 fármacos (pós-unificação case-insensitive)
    │  • group_by(Drug_Name_Norm)
    │  • mean_fc, median_fc, min_q_value, max(score)
    │  • direction = f(median_fc)
    │  • arrange(min_q, desc(score), desc(max_abs_fc))
    ▼
├── volcano_NQO1_drugs.png
├── dose_response_NQO1_drugs.png
├── ranking_drugs_NQO1.csv (completo)
├── top20_drugs_NQO1.csv
├── doxorubicin_NQO1_summary.csv
└── DGB_results_NQO1_from_xlsx_clean.csv
```

**STATUS:** ⚠️ Pipeline documentado corretamente, mas **o script NÃO foi re-executado** após as correções. Os outputs em disco ainda refletem o agrupamento case-sensitive antigo (1281 linhas, não 1266).

---

## ETAPA 2, AUDITORIA DO CÓDIGO (LINHA POR LINHA)

| # | Linha/Função | Problema | Gravidade |
|---|-------------|----------|-----------|
| 1 | `decimal_to_numeric` | Remove TODOS os pontos antes de converter vírgula. Se formato for "1.234.567,89" funciona, mas se for "1,234.56" (US), falha silenciosamente. Frágil. | BAIXA |
| 2 | `extract_dose_um` | Regex `\d+(?:\.\d+)?` espera ponto decimal. Depende de `str_replace_all(dose, ",", ".")` aplicado antes. OK, mas ordem frágil. | BAIXA |
| 3 | `read_nqo1_xlsx` | `expected_direction` calculado mas **nunca usado**. Código morto. | BAIXA |
| 4 | `rank_nqo1_drugs` | `score = max(score, na.rm = TRUE)`, usa o score MÁXIMO entre registros, não o score do efeito agregado. O paper descreve fórmula diferente. ⚠️ | ALTA |
| 5 | `rank_nqo1_drugs` | `direction = f(median_fc)`, OK, mas manumycin-a: mean_fc=-0,010, median_fc=+0,746. Contradição não sinalizada no código. | MÉDIA |
| 6 | `plot_dose_response` | Fallback quando `n_doses < 2` para todos os fármacos funciona, mas nunca foi testado (sempre há ≥1 fármaco com ≥2 doses). | BAIXA |
| 7 | `main()` | O filtro `doxo` usa `str_detect(Drug_Name_Norm, "doxo|...")`, OK, mas depende do script ser re-executado para ter `Drug_Name_Norm` no ranking. | MÉDIA |
| 8 | `main()` | Tabela 1 do paper (contagens por plataforma) **não é gerada pelo script**. O script não produz esta tabela. | ALTA |
| 9 | `cleanup_unneeded_outputs` | Lista hardcoded de arquivos. Se arquivos forem renomeados, a limpeza falha silenciosamente. | BAIXA |
| 10 | `AUDIT_EVENTS` | Cresce com `rbind()`, O(n²). Com ~15 eventos, irrelevante. | BAIXA |

---

## ETAPA 3, REPRODUÇÃO MENTAL DAS TABELAS

### Tabela 1 do manuscrito (Distribuição por plataforma)

| Claim no paper | Dado real (não-unificado, CSV atual) | Dado real (unificado, pós-script corrigido) | Status |
|---------------|--------------------------------------|---------------------------------------------|--------|
| LINCS: 1199 (94,7%) | 1199 entradas no ranking | 1206 fármacos únicos (95,3%) | 🔴 **ERRADO**, usa contagem não-unificada |
| CREEDS: 68 (5,4%) | 68 entradas | 68 fármacos (5,4%) | ✅ Coincide |
| CMAP: 7 (0,6%) | 7 entradas | 14 fármacos (1,1%) | 🔴 **ERRADO**, 7 são entradas não-unificadas |
| Multi-plataforma: 7 (0,6%) | 7 (CMAP+LINCS apenas) | **21** (14 CREEDS+LINCS, 6 CMAP+LINCS, 1 todos) | 🔴 **ERRADO**, 7 é apenas CMAP+LINCS |
| Down: 676 (53,4%) | 676 entradas não-unificadas | **668** (52,8%) | 🔴 **ERRADO**, 676+605=1281≠1266 |
| Up: 605 (47,8%) | 605 entradas não-unificadas | **598** (47,2%) | 🔴 **ERRADO**, mesmo problema |

**Conclusão:** A Tabela 1 inteira usa contagens do ranking NÃO-UNIFICADO (1281 entradas), mas as porcentagens são calculadas sobre o total UNIFICADO (1266). Isso é uma **inconsistência matemática**. A tabela precisa ser recalculada com dados unificados.

### Tabela 2 do manuscrito (Top 20)

| Fármaco | Paper (FC) | Paper (Score) | Unificado (FC) | Unificado (Score) | Status |
|---------|-----------|---------------|----------------|-------------------|--------|
| trichostatin A | +0,681 / 60,7 | 60,7 | +0,681 / 34,0 | 34,0 | 🔴 Score errado |
| Doxorubicina | −1,065 / 29,5 | 29,5 | −1,065 / 29,5 | 29,5 | ✅ (único unificado) |
| manumycin-a | −0,010 / 68,1 | 68,1 | −0,010 / 0,2 | 0,2 | 🔴 Score errado |
| parthenolide | +1,569 / 94,0 | 94,0 | +1,569 / 28,8 | 28,8 | 🔴 Score errado |
| panobinostat | −0,539 / 24,9 | 24,9 | −0,539 / 9,3 | 9,3 | 🔴 Score errado |
| piperlongumine | −0,705 / 49,0 | 49,0 | −0,705 / 12,0 | 12,0 | 🔴 Score errado |
| radicicol | +1,319 / 53,7 | 53,7 | +1,319 / 21,9 | 21,9 | 🔴 Score errado |
| BRD-K68548958 | −0,426 / 53,6 | 53,6 | −0,426 / 7,0 | 7,0 | 🔴 Score errado |
| YM-155 | −0,977 / 47,5 | 47,5 | −0,977 / 15,8 | 15,8 | 🔴 Score errado |
| COT-10b | −1,612 / 34,7 | 34,7 | −1,612 / 26,1 | 26,1 | 🔴 Score errado |
| 15-delta-PGJ2 | +1,468 / 34,8 | 34,8 | +1,468 / 23,6 | 23,6 | 🔴 Score errado |
| PHA-793887 | −0,941 / 57,4 | 57,4 | −0,941 / 15,0 | 15,0 | 🔴 Score errado |
| tacedinaline | +0,027 / 25,6 | 25,6 | +0,027 / 0,4 | 0,4 | 🔴 Score errado |
| SA-1447005 | +2,575 / 71,6 | 71,6 | +2,575 / 40,6 | 40,6 | 🔴 Score errado |
| CT-200783 | +1,261 / 67,0 | 67,0 | +1,261 / 19,6 | 19,6 | 🔴 Score errado |
| NSC-632839 | +0,682 / 53,5 | 53,5 | +0,682 / 10,6 | 10,6 | 🔴 Score errado |
| narciclasine | −1,539 / 52,3 | 52,3 | −1,539 / 23,6 | 23,6 | 🔴 Score errado |
| MW-A1-12 | −1,897 / 54,9 | 54,9 | −1,897 / 28,7 | 28,7 | 🔴 Score errado |
| withaferin-a | +0,552 / 43,1 | 43,1 | +0,552 / 8,2 | 8,2 | 🔴 Score errado |
| WZ-3105 | −1,623 / 51,4 | 51,4 | −1,623 / 23,6 | 23,6 | 🔴 Score errado |

**Conclusão:** 19 dos 20 scores na Tabela 2 estão ERRADOS. O paper usa `max(score)` do agrupamento não-unificado (que é `max(-log10(q) × |FC_individual|)`), não o score do efeito agregado. Apenas a doxorrubicina foi recalculada com a fórmula unificada. **A Tabela 2 é um híbrido inconsistente.**

### Implicações do erro de score:

- O paper descreve o score como `−log₁₀(q-value) × |Fold Change|` (fórmula do efeito agregado)
- Mas o código calcula `max(−log₁₀(q_i) × |FC_i|)` para cada registro individual i
- Isso significa que **a descrição metodológica não corresponde ao que o código faz**
- Os scores atuais refletem o MELHOR registro individual, não o efeito médio do fármaco
- Exemplo: trichostatin A tem 4 registros. O melhor tem |FC|=1,214 com q=0 → score=60,7. Mas o FC médio é 0,681 → score agregado = 34,0.

---

## ETAPA 4, AUDITORIA METODOLÓGICA

| Etapa declarada | Executada pelo script.R? | Evidência |
|----------------|--------------------------|-----------|
| Download GSE116250 | ❌ NÃO | Script externo não compartilhado |
| DEG com limma | ❌ NÃO | Script externo não compartilhado |
| GSEA/GSVA da via AMPK | ❌ NÃO | Script externo não compartilhado |
| Consulta DGB (geração do XLSX) | ❌ NÃO | Ferramenta web externa |
| Leitura e processamento do XLSX | ✅ SIM | `read_nqo1_xlsx()` |
| Ranking | ✅ SIM* | `rank_nqo1_drugs()`, mas não re-executado |
| Volcano plot | ✅ SIM | `plot_volcano()` |
| Dose-response | ✅ SIM | `plot_dose_response()` |
| Rede PPI (STRING) | ✅ SIM | `plot_ppi_network()` |
| PubChem | ⚠️ PARCIAL | Apenas carrega CSV pré-existente |
| Tabela 1 do paper | ❌ NÃO | **Não é gerada pelo script** |
| Tabela 2 do paper | ⚠️ PARCIAL | `top20_drugs_NQO1.csv` existe, mas com scores antigos |

**Falha crítica de reprodutibilidade:** A Tabela 1 (distribuição plataforma/direção) não é produzida por nenhuma função do script.R.

---

## ETAPA 5, AUDITORIA FRASE POR FRASE

| # | Frase | Classificação | Evidência |
|---|-------|--------------|-----------|
| 1 | "3523 registros válidos" | ✅ VERDADEIRA | audit_events: "Linhas validas: 3523" |
| 2 | "1266 fármacos únicos (pós-unificação)" | ✅ VERDADEIRA | 1266 unique lowercase names |
| 3 | "676 (53,4%) downreguladores e 605 (47,8%) upreguladores" | 🔴 FALSA | 676+605=1281≠1266. Unificado: 668 (52,8%) e 598 (47,2%) |
| 4 | "LINCS L1000 contribuiu com 94,7%... CREEDS 5,4%... CMAP 0,6%" | 🔴 FALSA | % calculadas sobre 1266 mas com numeradores não-unificados. Unificado: 95,3%, 5,4%, 1,1% |
| 5 | "Sete fármacos (0,6%) apresentaram evidência em mais de uma plataforma" | 🔴 FALSA | Unificado: 21 fármacos multi-plataforma (1,7%) |
| 6 | "trichostatin A (q-value < 10⁻⁵⁰; FC médio = +0,681; upregulador; CMAP)" | ✅ VERDADEIRA | CSV confirma |
| 7 | "parthenolide (q-value = 4,42×10⁻¹⁹; FC médio = +1,569; score = 94,0)" | ⚠️ SCORE ERRADO | q e FC corretos. Score unificado real = 28,8 |
| 8 | "doxorrubicina (q-value = 2,21×10⁻²⁸; FC médio = −1,065; downregulador; CREEDS e LINCS)" | ✅ VERDADEIRA | Dados unificados corretos |
| 9 | "n = 15, duas plataformas" | ✅ VERDADEIRA | 10 CREEDS + 5 LINCS = 15 |
| 10 | "11 (55%) downreguladores e 9 (45%) upreguladores" no top 20 | ⚠️ DEPENDE | Top 20 unificado pode ter proporção diferente |
| 11 | "19 dos 20 (95%) LINCS L1000" | ⚠️ DEPENDE | Unificado: 18/20 (90%), doxorubicina é CREEDS+LINCS |
| 12 | "escore composto definido como: score = −log₁₀(q-value) × |Fold Change|" | ⚠️ IMPRECISA | O código usa max(score_individual), não score_agregado |
| 13 | "O escore penaliza fármacos com q-values elevados" | ✅ VERDADEIRA | Conceitualmente correto |
| 14 | "direção consolidada foi determinada pelo Fold Change mediano" | ✅ VERDADEIRA | Confirmado no código |
| 15 | "trichostatin A... inibidor de histona desacetilases... ativação transcricional global" | ⚠️ SUPORTADA | Ref [8] suporta mecanismo de HDACi |
| 16 | "parthenolide... inibir NF-κB e induzir estresse oxidativo... via NRF2/ARE" | ⚠️ ESPECULATIVA | Ref [9] suporta NF-κB; conexão NRF2/ARE é inferência (ref [5] trata de NQO1/NRF2 mas não especificamente parthenolide) |
| 17 | "doxorrubicina reprime a expressão de NQO1... consistente com o perfil pró-oxidante" | ⚠️ PLAUSÍVEL | Observação + inferência mecanística. Ref [10] suporta perfil pró-oxidante |
| 18 | "A replicação deste achado em duas bases independentes confere maior confiabilidade" | ✅ VERDADEIRA | CREEDS + LINCS convergem em downregulation |
| 19 | "ranking... constitui um recurso sistemático e reproduzível" | ✅ VERDADEIRA | Código aberto, pipeline documentado |
| 20 | "Estudos futuros poderão validar experimentalmente os efeitos preditos" | ✅ APROPRIADA | Linguagem condicional correta |

---

## ETAPA 6, AUDITORIA FARMACOGENÔMICA

| Verificação | Resultado |
|-------------|-----------|
| Direção da doxorrubicina: DOWN | ✅ Consistente CREEDS + LINCS |
| Direção do trichostatin A: UP | ✅ CMAP (4 registros, todos up) |
| Direção do parthenolide: UP | ✅ LINCS (29 registros) |
| Direção do panobinostat: DOWN | ✅ LINCS (14 registros) |
| manumycin-a: UP por median_fc, mean_fc ≈ 0 | ⚠️ Documentado no paper, OK |
| Conflito entre bases? | Nenhum detectado para top 5 |
| Monensina | ✅ NÃO mencionada (corrigido) |
| Geldanamicina | ✅ NÃO mencionada como upreguladora (corrigido) |
| Wortmannin | ✅ NÃO mencionada como upreguladora (corrigido) |

---

## ETAPA 7, AUDITORIA ESTATÍSTICA

| Problema | Descrição | Gravidade |
|----------|-----------|-----------|
| **Score arbitrário** | `score = max(−log10(q) × |FC|)`. Sem validação, sem justificativa estatística formal, sem comparação com alternativas (ex: weighted mean, meta-analysis). | ALTA |
| **Divergência score descrito vs. implementado** | O paper diz `score = −log10(q) × |FC|` (agregado). O código faz `max(−log10(q_i) × |FC_i|)` (melhor registro). São conceitos diferentes. | ALTA |
| **Dependência entre registros** | Múltiplos registros do mesmo fármaco na mesma plataforma (diferentes doses/tempos) NÃO são independentes. A significância pode estar inflada. | MÉDIA |
| **Desbalanceamento de plataformas** | 94,7% dos registros são LINCS L1000 (~1000 genes). CMAP (microarranjo completo, ~22000 genes) tem apenas 7 fármacos. Viés de representação. | MÉDIA |
| **q-value = 0** | trichostatin A: q=0 capado em 10⁻⁵⁰. O valor real de -log10(q) é ∞. O cap é arbitrário. | BAIXA |
| **Sem correção para múltiplos testes** | 1266 fármacos testados. Os q-values do DGB já incorporam correção? Não está documentado. | MÉDIA |
| **n_tests variável** | Fármacos com 1 registro (n=1) vs. 29 registros (parthenolide) recebem o mesmo peso no ranking. | MÉDIA |

---

## ETAPA 8, AUDITORIA BIOLÓGICA

| Interpretação no paper | Problema |
|------------------------|----------|
| "ativação transcricional global [por HDACi], o que é consistente com a upregulação de NQO1" | ✅ Correlação plausível com suporte da literatura |
| "poderia explicar a indução de NQO1 como resposta compensatória celular mediada pela via NRF2/ARE" | ⚠️ **ESPECULATIVA.** Mecanismo NRF2/ARE → NQO1 existe, mas NÃO há dado neste estudo que demonstre que o parthenolide induz NQO1 VIA NRF2. É inferência. |
| "doxorrubicina reprime... consistente com o perfil pró-oxidante" | ⚠️ **PLAUSÍVEL mas não demonstrado.** Correlação entre "droga pró-oxidante" e "reprime gene antioxidante" é observacional. |
| "A replicação deste achado em duas bases independentes confere maior confiabilidade" | ✅ CORRETO. Duas plataformas independentes com mesma direção. |

**Todas as interpretações mecanísticas permanecem como HIPÓTESES.** Nenhuma foi validada experimentalmente. O paper reconhece isso nas limitações.

---

## ETAPA 9, AUDITORIA DAS REFERÊNCIAS

| Ref | Suporta? | Observação |
|-----|----------|------------|
| [1] HODOS 2016 | ✅ | In silico drug repurposing |
| [2] EKINS 2019 | ⚠️ | ML para drug discovery, tangencial |
| [3] LAMB 2006 | ✅ | Connectivity Map original |
| [4] ROSS 2017 | ✅ | NQO1 functions |
| [5] DINKOVA-KOSTOVA 2010 | ✅ | NQO1/NRF2 cytoprotection |
| [6] WANG 2021 | ✅ | DGB platform |
| [7] SUBRAMANIAN 2017 | ✅ | LINCS L1000 platform |
| [8] YOSHIDA 1995 | ⚠️ | Trichostatin A/HDAC, **1995, desatualizada**. Existem revisões mais recentes. |
| [9] GHAZNAVI 2023 | ✅ | Parthenolide review |
| [10] WALLACE 2003 | ⚠️ | Doxorubicin cardiotoxicity, **2003, 23 anos**. Existem revisões mais recentes. |
| [11] BOZKURT 2021 | ⚠️ | CMD review, citada apenas como contexto. OK. |
| [12] ZANGIROLAMI-RAIMUNDO 2018 | ⚠️ | Estudos transversais, **tangencial**. Não suporta nenhuma afirmação específica. |

**Recomendações:**
- Ref [8]: Substituir por revisão mais recente de HDACi (ex: Seto & Yoshida, Cold Spring Harb Perspect Med, 2014)
- Ref [10]: Substituir por revisão mais recente de doxorubicina (ex: Rawat et al., 2021 ou Christidi & Brunham, 2023)
- Ref [12]: Remover ou justificar melhor

---

## ETAPA 10, AUDITORIA EDITORIAL

| Aspecto | Avaliação |
|---------|-----------|
| Clareza | ✅ Boa. Texto claro e objetivo. |
| Coesão | ✅ Boa. Seções fluem logicamente. |
| Coerência | ⚠️ Dados da Tabela 1 inconsistentes com a narrativa. |
| Repetições | ⚠️ "A doxorrubicina é destacada com destaque visual", redundante. |
| Jargão excessivo | ✅ OK para público-alvo. |
| Linguagem causal | ✅ Controlada. Uso de "consistente com", "poderia explicar", "sugere". |
| Linguagem translacional | ✅ Removida (corrigido da versão anterior). |
| Linguagem clínica | ✅ Removida (corrigido da versão anterior). |

---

## ETAPA 11, PARECER COMO REVISOR DA NATURE

### Major Comments

**1. Inconsistência crítica nos dados da Tabela 1.**
As contagens de fármacos por direção (676 + 605 = 1281) não correspondem ao total declarado de 1266 fármacos únicos. As porcentagens (53,4% e 47,8%) são calculadas com denominador unificado (1266) mas numeradores não-unificados (1281). Além disso, o número de fármacos multi-plataforma (7) está subestimado, a unificação case-insensitive revela 21 fármacos com evidência em múltiplas bases. A Tabela 1 deve ser completamente recalculada após a execução do script corrigido com agrupamento case-insensitive.

**2. Divergência entre o escore descrito e o escore implementado.**
O manuscrito define o escore como `−log₁₀(q-value) × |Fold Change|`, mas o código (script.R, função `rank_nqo1_drugs`) calcula `max(−log₁₀(qᵢ) × |FCᵢ|)` para registros individuais, não para o efeito agregado. Isso produz scores substancialmente diferentes (ex: trichostatin A: 60,7 pelo código vs. 34,0 pela fórmula descrita). Os autores devem alinhar a descrição metodológica com a implementação real ou alterar o código.

**3. A Tabela 2 é um híbrido metodologicamente inconsistente.**
Dezenove dos 20 fármacos na Tabela 2 usam scores do agrupamento não-unificado (case-sensitive). Apenas a doxorrubicina foi recalculada com valores unificados. Isso torna a tabela internamente inconsistente e potencialmente enganosa.

**4. A Tabela 1 não é gerada pelo script fornecido.**
O script.R não produz a tabela de distribuição plataforma/direção apresentada como Tabela 1, comprometendo a reprodutibilidade.

### Minor Comments

**5.** Referências [8] (Yoshida 1995) e [10] (Wallace 2003) estão desatualizadas. Existem revisões mais recentes e abrangentes.

**6.** A frase "A doxorrubicina é destacada com destaque visual" contém redundância.

**7.** A referência [12] (Zangirolami-Raimundo 2018, estudos transversais) é tangencial e não suporta diretamente nenhuma afirmação do manuscrito.

**8.** O tratamento de q-value = 0 (capado em 10⁻⁵⁰) é arbitrário. Os autores devem justificar a escolha deste valor de corte.

### Fatal Flaws

- **F1:** Tabela 1 com inconsistência matemática (676+605≠1266)
- **F2:** Divergência entre escore descrito e escore implementado
- **F3:** Tabela 2 híbrida (19 scores não-unificados + 1 unificado)

### Mandatory Revisions

1. Re-executar o script.R corrigido e atualizar TODOS os outputs
2. Recalcular Tabela 1 com contagens unificadas
3. Recalcular Tabela 2 com scores consistentes (todos unificados)
4. Alinhar descrição metodológica do escore com a implementação
5. Gerar Tabela 1 programaticamente no script.R
6. Atualizar referências [8] e [10]

---

## ETAPA 12, TESTE DE ROBUSTEZ

> *"Se eu remover esta conclusão, o estudo continua válido?"*

| Conclusão | Remove? | Justificativa |
|-----------|---------|---------------|
| "trichostatin A, doxorrubicina, manumycin-a, parthenolide e panobinostat são os 5 com maior evidência estatística" | ✅ MANTÉM | Diretamente dos dados |
| "SA-1447005, MW-A1-12 e parthenolide têm maior magnitude" | ✅ MANTÉM | Diretamente dos dados |
| "doxorrubicina tem evidência convergente em 2 plataformas" | ✅ MANTÉM | Diretamente dos dados |
| "ranking pode orientar validação experimental futura" | ✅ MANTÉM | Conclusão apropriada para estudo in silico |
| "ativação transcricional global... consistente com upregulação" | ⚠️ QUESTIONÁVEL | Inferência mecanística. Remover não afeta validade do ranking. |
| "via NRF2/ARE" para parthenolide | ⚠️ REMOVER | Especulativa. Sem dados sobre NRF2/ARE neste estudo. |
| "consistente com perfil pró-oxidante" da doxorrubicina | ⚠️ QUESTIONÁVEL | Plausível mas inferencial. Pode ser mantida como hipótese. |

---

## ETAPA 13, 🔒 LOCK FINAL

| Item | Definição |
|------|-----------|
| **Objetivo definitivo** | Identificar, caracterizar e ranquear fármacos com evidência transcriptômica de modulação da expressão do gene NQO1, integrando as bases CMAP, LINCS L1000 e CREEDS via DGB. |
| **Pergunta científica única** | Quais fármacos apresentam evidência transcriptômica de modular a expressão de NQO1, e qual sua posição relativa no ranking integrado? |
| **Hipótese** | NULL (estudo descritivo-exploratório). Não há hipótese a ser testada; o objetivo é gerar um ranking reproduzível. |
| **Escopo** | Exclusivamente in silico. Triagem farmacogenômica de fármacos moduladores de NQO1 usando dados públicos (CMAP, LINCS, CREEDS). |
| **Dados utilizados** | DGB_results_NQO1.xlsx (6 abas); ranking_drugs_NQO1.csv; DGB_results_NQO1_from_xlsx_clean.csv |
| **Dados NÃO utilizados** | GSE116250 (apenas contexto); STRING (apenas PPI contextual); PubChem (output vazio); DEG_full_table.csv (apenas para identificar NQO1) |
| **Métodos utilizados** | Padronização de nomenclatura (case-insensitive); conversão de unidades (→µM); escore composto; agrupamento por fármaco; ranking por q-value e score |
| **Métodos NÃO utilizados** | GEOquery; limma; GSEA; GSVA; consulta PubChem ativa; qualquer validação experimental |
| **Limitações obrigatórias** | (1) In silico, requer validação experimental; (2) Dados de linhagens tumorais, podem não refletir tipos celulares primários; (3) Predominância de LINCS L1000 (genes Landmark, ~1000 transcritos); (4) Escore composto não validado externamente; (5) Variabilidade experimental entre registros do mesmo fármaco |
| **Conclusões PERMITIDAS** | (a) Ranking dos 5 fármacos com maior evidência estatística; (b) Fármacos com maior magnitude de modulação; (c) Distribuição up/down; (d) Distribuição por plataforma; (e) Fármacos com evidência multi-plataforma; (f) Recurso para priorização de validação experimental futura |
| **Conclusões PROIBIDAS** | ❌ Mecanismos moleculares; ❌ Causalidade; ❌ Implicações clínicas; ❌ Farmacovigilância; ❌ Enfermagem de precisão; ❌ Cardiotoxicidade demonstrada; ❌ Validação da via AMPK; ❌ Proteção redox; ❌ Resposta adaptativa/compensatória |
| **Título definitivo** | **Ranqueamento Farmacogenômico in silico de Fármacos Moduladores da Expressão do Gene NQO1 por Integração das Bases CMAP, LINCS L1000 e CREEDS** |
| **Estrutura definitiva** | Resumo (objetivo→método→resultados→conclusão) → Introdução (farmacogenômica computacional, NQO1, DGB, objetivo) → Métodos (origem, padronização, escore, ranking, código) → Resultados (Tabela 1 distribuição, Tabela 2 top 20, visualização) → Discussão (interpretação dos top fármacos, limitações) → Conclusão (ranking, magnitudes, recurso futuro) |

---

## 📊 RESUMO DE GRAVIDADE

| Gravidade | # Problemas | Itens |
|-----------|------------|-------|
| 🔴 **CRÍTICA** | 3 | Tabela 1 inconsistente (676+605≠1266); Escore descrito ≠ implementado; Tabela 2 híbrida |
| 🟠 **ALTA** | 4 | Contagens multi-plataforma erradas; Tabela 1 não gerada pelo script; Scores não-unificados; Referências desatualizadas |
| 🟡 **MÉDIA** | 5 | Dependência entre registros; Desbalanceamento de plataformas; Correção múltipla não documentada; n_tests variável ignorado; Interpretações mecanísticas especulativas |
| 🟢 **BAIXA** | 4 | Redundância textual; Código morto (expected_direction); Referência tangencial [12]; Cap arbitrário em q=0 |

---

## ✅ AÇÕES CORRETIVAS OBRIGATÓRIAS

1. **Re-executar `script.R`** com as correções de case-insensitivity → gerar outputs unificados
2. **Recalcular Tabela 1** com dados unificados (668 down, 598 up; 1206 LINCS, 68 CREEDS, 14 CMAP; 21 multi-plataforma)
3. **Recalcular Tabela 2** inteira com scores unificados consistentes
4. **Alinhar seção 2.3** (descrição do escore) com a implementação real OU alterar o código
5. **Adicionar geração da Tabela 1 ao script.R**
6. **Atualizar referências** [8] e [10] para versões mais recentes
7. **Remover** referência [12] ou justificar
8. **Corrigir** "destacada com destaque visual"
9. **Qualificar** interpretações mecanísticas como hipóteses não testadas

---

*Relatório produzido em 2026-06-26. Auditoria forense definitiva, LOCK FINAL.*
