#!/usr/bin/env Rscript

# ==============================================================================
# 05 — GSEA E GSVA DA VIA AMPK (versão didática, passo a passo)
# ==============================================================================
#
# ESTE SCRIPT É OPCIONAL.
#   Ele continua a análise depois do volcano. Se a aula termina no volcano,
#   você não precisa rodá-lo.
#
# O QUE É GSEA?
#   Gene Set Enrichment Analysis. Em vez de testar gene por gene com um corte
#   de p-valor, o GSEA ordena TODOS os genes por log2FC e pergunta:
#   "os genes da via AMPK estão concentrados no topo ou na base desse ranking?".
#
# O QUE É GSVA?
#   Gene Set Variation Analysis. O GSEA dá um resultado para o GRUPO. A GSVA
#   transforma a via em UM ESCORE POR AMOSTRA, como se a via fosse um
#   "supergene". Assim podemos comparar DCM vs CTRL amostra a amostra.
#
# COMO ESTE SCRIPT ESTÁ ORGANIZADO?
#   PARTE 1 — carregar dados e a via AMPK
#   PARTE 2 — GSEA (ranking -> teste -> tabela -> gráfico)
#   PARTE 3 — GSVA (escore por amostra -> teste -> tabela -> gráfico)
# ==============================================================================

# ------------------------------------------------------------------------------
# 0) CARREGAR A BASE (configuração, caminhos e funções auxiliares)
# ------------------------------------------------------------------------------
pasta_do_script <- function() {
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(args)) {
    return(dirname(normalizePath(sub("^--file=", "", args[1]),
                                  winslash = "/", mustWork = FALSE)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(p)) {
      return(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    }
  }
  getwd()
}
SCRIPT_DIR <- pasta_do_script()
if (!exists("preparar_dados", mode = "function")) {
  source(file.path(SCRIPT_DIR, "00_configuracao_e_dados.R"), encoding = "UTF-8")
}

suppressPackageStartupMessages({
  library(fgsea)     # GSEA
  library(GSVA)      # GSVA
  library(dplyr)     # manipulação de tabelas
  library(readr)     # leitura/escrita de CSV
  library(ggplot2)   # gráficos
})

# ==============================================================================
# PARTE 1 — CARREGAR OS DADOS E A VIA AMPK
# ==============================================================================
step("Carregando dados e a via AMPK")

# O resultado do limma foi criado pelo script 02_volcano.R.
arquivo_limma <- file.path(PROCESSED, "resultado_limma_DCM_vs_CTRL.rds")
if (!file.exists(arquivo_limma)) {
  stop("Não encontrei o resultado do limma.\n",
       "Rode primeiro o script 02_volcano.R.", call. = FALSE)
}

resultado_limma <- readRDS(arquivo_limma)

# Separamos os três objetos que vamos usar:
#   deg      = tabela com gene, log2FoldChange, padj etc.
#   expr_de  = matriz de expressão log2 (genes x amostras)
#   pheno_de = tabela de grupos das amostras (CTRL/DCM)
deg      <- resultado_limma$deg
expr_de  <- resultado_limma$expr_de
pheno_de <- resultado_limma$pheno_de

# Lista de genes da via AMPK (usa cache; não baixa de novo se já existe).
ampk <- baixar_genes_ampk()

# Recorte da via dentro da tabela de expressão diferencial.
deg_ampk <- deg |>
  filter(gene %in% ampk$symbol)

message("Genes da via AMPK no ranking: ", nrow(deg_ampk))

# ==============================================================================
# PARTE 2 — GSEA
# ==============================================================================
step("GSEA da via AMPK")

# ------------------------------------------------------------------------------
# 2.1) MONTAR O RANKING
# ------------------------------------------------------------------------------
# O GSEA precisa de uma lista ORDENADA de genes. Ordenamos do maior log2FC para
# o menor:
#   topo = genes mais "up" na DCM
#   base = genes mais "down" na DCM
#
# NÃO aplicamos corte de p-valor. TODOS os genes entram no ranking. É justamente
# isso que diferencia o GSEA de uma lista de genes significativos.
ranking <- deg |>
  filter(!is.na(log2FoldChange)) |>
  arrange(desc(log2FoldChange)) |>
  mutate(posicao = row_number()) |>
  select(posicao, gene, log2FoldChange)

# O fgsea pede um vetor numérico com NOMES. Aqui:
#   valores = log2FoldChange
#   nomes   = símbolo do gene
ranking_log2fc <- ranking$log2FoldChange
names(ranking_log2fc) <- ranking$gene

message("Genes no ranking: ", length(ranking_log2fc))
message("Topo do ranking: ",
        paste(head(ranking$gene, 5), collapse = ", "))
message("Base do ranking: ",
        paste(tail(ranking$gene, 5), collapse = ", "))

# ------------------------------------------------------------------------------
# 2.2) RODAR O GSEA
# ------------------------------------------------------------------------------
# fgseaSimple() testa se os genes da via AMPK se concentram nas pontas do
# ranking. nperm = 10.000 permutações são usadas para calcular a significância.
#
# O resultado tem uma linha por via testada. Como testamos só a via AMPK,
# esperamos uma linha.
#
# set.seed() imediatamente antes do fgseaSimple() garante que o resultado das
# permutações não dependa do estado do gerador aleatório deixado pelos scripts
# anteriores (por exemplo, o layout da rede no script 04).
set.seed(42)
resultado_fgsea <- fgseaSimple(
  pathways = list(AMPK = ampk$symbol),
  stats    = ranking_log2fc,
  nperm    = 10000
)

if (nrow(resultado_fgsea) == 0) {
  stop("O GSEA não retornou resultado. Verifique se os símbolos da via AMPK ",
       "estão presentes no ranking.", call. = FALSE)
}

# Pegamos a primeira (e única) linha de forma explícita, para ficar claro.
resultado_gsea <- resultado_fgsea[1, ]

# ------------------------------------------------------------------------------
# 2.3) ORGANIZAR A TABELA DE RESULTADO
# ------------------------------------------------------------------------------
# ES  = enrichment score. >0 = via no topo (up); <0 = via na base (down).
# NES = ES normalizado pelo tamanho da via. Permite comparar vias.
# p_value / p_ajustado = significância.
# leadingEdge = genes da via que mais contribuíram para o sinal.
tabela_gsea <- data.frame(
  via               = "KEGG_AMPK_SIGNALING_PATHWAY_HSA04152",
  n_genes           = resultado_gsea$size,
  ES                = resultado_gsea$ES,
  NES               = resultado_gsea$NES,
  p_value           = resultado_gsea$pval,
  p_ajustado        = resultado_gsea$padj,
  genes_leading_edge = paste(resultado_gsea$leadingEdge[[1]], collapse = "; "),
  stringsAsFactors  = FALSE
)
write_csv(tabela_gsea, file.path(TABLES, "gsea_via_ampk.csv"))

# Tabela só com os genes do leading edge e o log2FC de cada um. Ajuda a
# entender QUAIS genes puxaram o sinal da via.
genes_leading_edge <- resultado_gsea$leadingEdge[[1]]
tabela_leading_edge <- data.frame(
  gene = genes_leading_edge,
  log2FoldChange = ranking$log2FoldChange[match(genes_leading_edge, ranking$gene)],
  stringsAsFactors = FALSE
)
write_csv(tabela_leading_edge,
          file.path(TABLES, "gsea_via_ampk_leading_edge.csv"))

message(sprintf("GSEA AMPK: NES = %.3f | p = %.4f | tamanho = %d genes",
                resultado_gsea$NES, resultado_gsea$pval, resultado_gsea$size))
message("Leading edge (primeiros 10): ",
        paste(head(genes_leading_edge, 10), collapse = ", "))

# ------------------------------------------------------------------------------
# 2.4) GRÁFICO DE ENRIQUECIMENTO
# ------------------------------------------------------------------------------
# COMO LER A CURVA?
#   * O eixo x é o ranking dos genes (esquerda = down, direita = up na DCM).
#   * Cada traço vertical é UM gene da via AMPK.
#   * A linha sobe quando encontramos genes da via e desce quando passamos por
#     genes que não são da via. O ponto mais alto (ou mais baixo) é o ES.
#   * Curva subindo no lado direito = via enriquecida em genes up.
#   * Curva descendo no lado esquerdo = via enriquecida em genes down.
grafico_gsea <- plotEnrichment(ampk$symbol, ranking_log2fc) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    title    = "GSEA da via AMPK (hsa04152)",
    subtitle = sprintf("NES = %.2f | p = %.4f | tamanho da via = %d genes",
                       resultado_gsea$NES, resultado_gsea$pval,
                       resultado_gsea$size),
    x = "Ranking dos genes por log2FC (down na DCM <- -> up na DCM)",
    y = "Enrichment score (ES)",
    caption = paste0("Cada traço vertical é um gene da via. ",
                     "ES > 0 = via concentrada no topo (up); ",
                     "ES < 0 = na base (down).")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.caption = element_text(hjust = 0))

save_plot(grafico_gsea, "GSEA_AMPK_enrichment_plot.png", 8, 6)

# ==============================================================================
# PARTE 3 — GSVA
# ==============================================================================
step("GSVA da via AMPK")

# ------------------------------------------------------------------------------
# 3.1) PREPARAR A MATRIZ E A LISTA DE GENES
# ------------------------------------------------------------------------------
# A GSVA espera:
#   * uma matriz de expressão (genes nas linhas, amostras nas colunas);
#   * uma lista de conjuntos de genes (aqui, só a via AMPK).
matriz_expressao <- as.matrix(expr_de)
conjuntos_genes  <- list(AMPK = ampk$symbol)

# kcdf = "Gaussian" porque os dados já são contínuos em log2 (RNA-seq).
parametros_gsva <- gsvaParam(
  exprData   = matriz_expressao,
  geneSets   = conjuntos_genes,
  kcdf       = "Gaussian"
)

# ------------------------------------------------------------------------------
# 3.2) RODAR A GSVA
# ------------------------------------------------------------------------------
# O resultado é uma matriz 1 x n_amostras: um escore da via para cada amostra.
matriz_escores <- gsva(parametros_gsva)

# ------------------------------------------------------------------------------
# 3.3) MONTAR A TABELA POR AMOSTRA
# ------------------------------------------------------------------------------
# Cada linha = uma amostra; cada coluna = nome da amostra, escore e grupo.
tabela_escores <- data.frame(
  sample     = colnames(matriz_escores),
  AMPK_score = as.numeric(matriz_escores["AMPK", ]),
  group      = pheno_de$group[match(colnames(matriz_escores), pheno_de$sample)],
  stringsAsFactors = FALSE
)
write_csv(tabela_escores, file.path(TABLES, "gsva_via_ampk.csv"))

# ------------------------------------------------------------------------------
# 3.4) RESUMO POR GRUPO
# ------------------------------------------------------------------------------
# Antes do teste, vale olhar as estatísticas descritivas. Isso ajuda a entender
# se a diferença entre os grupos é grande ou pequena.
resumo_grupos <- tabela_escores |>
  group_by(group) |>
  summarise(
    n       = n(),
    media   = mean(AMPK_score),
    mediana = median(AMPK_score),
    desvio  = sd(AMPK_score),
    .groups = "drop"
  )
write_csv(resumo_grupos, file.path(TABLES, "gsva_via_ampk_resumo_grupos.csv"))

cat("\nResumo do escore GSVA por grupo:\n")
print(as.data.frame(resumo_grupos))

# ------------------------------------------------------------------------------
# 3.5) TESTE ESTATÍSTICO
# ------------------------------------------------------------------------------
# Wilcoxon é um teste não paramétrico: compara a distribuição dos escores entre
# DCM e CTRL sem assumir normalidade. É adequado para amostras pequenas.
teste_wilcox <- wilcox.test(AMPK_score ~ group, data = tabela_escores)

p_wilcox     <- teste_wilcox$p.value
p_wilcox_log <- -log10(p_wilcox)   # escala log: quanto maior, mais significativo

message(sprintf("GSVA Wilcoxon DCM vs CTRL: p = %.4g | -log10(p) = %.2f",
                p_wilcox, p_wilcox_log))

# Tabela-resumo do teste. Além do p-valor bruto, guardamos o -log10(p), que é
# a forma mais usada em gráficos de omics (valores maiores = mais significativos).
tabela_teste_gsva <- data.frame(
  comparacao  = "DCM vs CTRL",
  teste       = "Wilcoxon",
  p_value     = p_wilcox,
  neg_log10_p = p_wilcox_log,
  n_CTRL      = sum(tabela_escores$group == "CTRL"),
  n_DCM       = sum(tabela_escores$group == "DCM"),
  stringsAsFactors = FALSE
)
write_csv(tabela_teste_gsva,
          file.path(TABLES, "gsva_via_ampk_teste.csv"))

# ------------------------------------------------------------------------------
# 3.6) GRÁFICO
# ------------------------------------------------------------------------------
# COMO LER?
#   * Cada ponto é uma amostra.
#   * A caixa mostra mediana e quartis.
#   * Se as caixas de CTRL e DCM estão em alturas diferentes, a via tem
#     atividade média diferente entre os grupos.
#   * O p-valor e o -log10(p) aparecem no subtítulo.
grafico_gsva <- ggplot(tabela_escores, aes(x = group, y = AMPK_score, fill = group)) +
  geom_boxplot(width = 0.55, alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.7, color = "grey25") +
  scale_fill_manual(values = c(CTRL = "#1B9E77", DCM = "#D95F02")) +
  labs(
    title    = "GSVA da via AMPK (escore por amostra)",
    subtitle = sprintf("Wilcoxon: p = %.3g | -log10(p) = %.2f | n = %d amostras",
                       p_wilcox, p_wilcox_log, nrow(tabela_escores)),
    x = NULL,
    y = "Escore GSVA da via AMPK",
    fill = "Grupo",
    caption = paste0("Cada ponto é uma amostra; a caixa mostra a mediana e os quartis. ",
                     "Escore alto = via mais ativa naquela amostra.")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0))

save_plot(grafico_gsva, "GSVA_AMPK_boxplot.png", 7.5, 6)

# ==============================================================================
# FIM
# ==============================================================================
step("FIM DO SCRIPT 05 — GSEA e GSVA concluídos.")
message("GSEA: output/figures/GSEA_AMPK_enrichment_plot.png")
message("GSVA: output/figures/GSVA_AMPK_boxplot.png")
message("Tabelas: output/tables/gsea_via_ampk*.csv e gsva_via_ampk*.csv")
