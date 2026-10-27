#!/usr/bin/env Rscript

# ==============================================================================
# 06 — ROBUSTEZ: CORREÇÃO DE BATCH (SVA) E DECONVOLUÇÃO CELULAR (xCell)
# ==============================================================================
#
# O QUE ESTE SCRIPT FAZ?
#   1) Aplica uma correção de efeitos latentes/batch com SVA e compara o
#      resultado do limma com e sem a correção.
#   2) Estima a composição celular aproximada do tecido com xCell.
#
# POR QUE SVA?
#   O GSE116250 não traz uma coluna explícita de lote no series matrix.
#   O SVA (Surrogate Variable Analysis) estima fatores latentes a partir da
#   própria matriz de expressão, sem usar os rótulos de grupo para encontrá-los.
#   Incluímos dois fatores no modelo para testar se o efeito do NQO1 é robusto.
#
# POR QUE xCell?
#   O tecido cardíaco é uma mistura de cardiomiócitos, fibroblastos, células
#   endoteliais, células imunes e outros tipos. O xCell converte a expressão
#   bulk em escores de enriquecimento por tipo celular. Os escores são
#   relativos e não devem ser lidos como frações absolutas.
#
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(42)

# ------------------------------------------------------------------------------
# 0) CARREGAR A BASE
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
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(limma)
  library(sva)
  library(xCell)
})

# ------------------------------------------------------------------------------
# 1) DADOS E RESULTADO DO LIMMA SEM CORREÇÃO
# ------------------------------------------------------------------------------
dados <- preparar_dados()
expr  <- dados$expr
pheno <- dados$pheno

rds_limma <- file.path(PROCESSED, "resultado_limma_DCM_vs_CTRL.rds")
if (!file.exists(rds_limma)) {
  stop("Rode primeiro o script 02_volcano.R para gerar o resultado do limma.",
       call. = FALSE)
}
resultado <- readRDS(rds_limma)
deg       <- resultado$deg
pheno_de  <- resultado$pheno_de
expr_de   <- resultado$expr_de

# ------------------------------------------------------------------------------
# 2) CORREÇÃO DE BATCH/EFEITOS LATENTES COM SVA
# ------------------------------------------------------------------------------
# O design principal preserva grupo, sexo e idade. O SVA procura fatores que
# explicam variabilidade não atribuída a essas variáveis. Usamos dois fatores
# para uma correção conservadora (o critério automático pode sugerir muitos).
mod  <- model.matrix(~ 0 + group + sex + age, data = pheno_de)
colnames(mod) <- sub("^group", "", colnames(mod))
mod0 <- model.matrix(~ sex + age, data = pheno_de)

set.seed(42)
sv <- sva::sva(expr_de, mod, mod0, n.sv = 2)
message("SVA: ", sv$n.sv, " fatores latentes incluídos no modelo.")

design_sv <- cbind(mod, sv$sv)
colnames(design_sv)[(ncol(mod) + 1):ncol(design_sv)] <-
  paste0("SV", seq_len(sv$n.sv))

contraste_sv <- makeContrasts("DCM-CTRL", levels = design_sv)
fit_sv <- eBayes(
  contrasts.fit(lmFit(expr_de, design_sv), contraste_sv),
  trend = TRUE
)
tt_sv <- topTable(fit_sv, number = Inf, sort.by = "none")

deg_sv <- data.frame(
  gene           = rownames(tt_sv),
  log2FoldChange = tt_sv$logFC,
  AveExpr        = tt_sv$AveExpr,
  pvalue         = tt_sv$P.Value,
  padj           = tt_sv$adj.P.Val,
  stringsAsFactors = FALSE
) |>
  mutate(
    direction   = ifelse(log2FoldChange > 0, "up", "down"),
    significant = padj < 0.05
  )

write_csv(deg_sv, file.path(TABLES, "expressao_diferencial_com_SVA.csv"))

fatores <- data.frame(
  sample = pheno_de$sample,
  group  = pheno_de$group,
  sv$sv,
  check.names = FALSE
)
colnames(fatores)[3:ncol(fatores)] <- paste0("SV", seq_len(sv$n.sv))
write_csv(fatores, file.path(TABLES, "SVA_fatores.csv"))

# Comparação antes/depois
nq_sem <- deg[deg$gene == "NQO1", ]
nq_com <- deg_sv[deg_sv$gene == "NQO1", ]

comparacao <- data.frame(
  metrica = c(
    "NQO1 log2FC",
    "NQO1 padj",
    "DEGs padj < 0,05",
    "DEGs padj < 0,05 e |log2FC| >= 0,585",
    "Correlação log2FC (todos os genes)"
  ),
  sem_SVA = c(
    nq_sem$log2FoldChange,
    nq_sem$padj,
    sum(deg$padj < 0.05, na.rm = TRUE),
    sum(deg$padj < 0.05 & abs(deg$log2FoldChange) >= 0.585, na.rm = TRUE),
    1
  ),
  com_SVA = c(
    nq_com$log2FoldChange,
    nq_com$padj,
    sum(deg_sv$padj < 0.05, na.rm = TRUE),
    sum(deg_sv$padj < 0.05 & abs(deg_sv$log2FoldChange) >= 0.585, na.rm = TRUE),
    cor(deg$log2FoldChange[match(deg_sv$gene, deg$gene)],
        deg_sv$log2FoldChange, use = "complete.obs")
  )
)
write_csv(comparacao, file.path(TABLES, "SVA_comparacao.csv"))

# Gráfico antes/depois
sc <- data.frame(
  gene = deg_sv$gene,
  sem  = deg$log2FoldChange[match(deg_sv$gene, deg$gene)],
  com  = deg_sv$log2FoldChange
)

p_sva <- ggplot(sc, aes(sem, com)) +
  geom_point(alpha = 0.25, size = 0.9, color = "grey45") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey60") +
  geom_point(data = filter(sc, gene == "NQO1"),
             color = "#D95F02", size = 3.5) +
  geom_text_repel(data = filter(sc, gene == "NQO1"),
                  aes(label = gene), fontface = "bold",
                  color = "#D95F02", size = 5, seed = 42) +
  labs(
    title = "Efeito da correção de batch (SVA) no log2FC",
    subtitle = paste0("Cada ponto é um gene; a linha tracejada é y = x. ",
                      "NQO1 está destacado."),
    x = "log2FC sem SVA",
    y = "log2FC com SVA"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0))

save_plot(p_sva, "SVA_logFC_antes_depois.png", 7.5, 6.2)

# ------------------------------------------------------------------------------
# 3) DECONVOLUÇÃO CELULAR COM xCell
# ------------------------------------------------------------------------------
# Usamos todas as 64 amostras para o mapa geral; a comparação estatística é
# feita apenas em CTRL vs DCM, como no limma.
expr_log <- log2(expr + 1)
expr_log <- expr_log[rowSums(expr > 1) >= 0.25 * ncol(expr), , drop = FALSE]

message("Rodando xCell (pode levar alguns minutos)...")
set.seed(42)
xc <- xCell::xCellAnalysis(expr_log, rnaseq = TRUE, parallel.sz = 1)

xc_df <- as.data.frame(xc) |>
  tibble::rownames_to_column("cell_type")
write_csv(xc_df, file.path(TABLES, "xCell_scores.csv"))

# Teste DCM vs CTRL para cada tipo celular
xc_de <- xc[, pheno_de$sample, drop = FALSE]
teste_xcell <- data.frame(
  cell_type = rownames(xc_de),
  CTRL      = rowMeans(xc_de[, pheno_de$group == "CTRL", drop = FALSE]),
  DCM       = rowMeans(xc_de[, pheno_de$group == "DCM", drop = FALSE]),
  stringsAsFactors = FALSE
) |>
  mutate(diff = DCM - CTRL)

teste_xcell$p <- apply(xc_de, 1, function(v) {
  wilcox.test(v[pheno_de$group == "CTRL"],
              v[pheno_de$group == "DCM"])$p.value
})
teste_xcell$padj <- p.adjust(teste_xcell$p, method = "BH")
teste_xcell <- teste_xcell |>
  arrange(padj, desc(abs(diff)))

write_csv(teste_xcell, file.path(TABLES, "xCell_teste.csv"))

# Tipos celulares relevantes para o coração e o microambiente
celulas <- c(
  "Myocytes", "Fibroblasts", "Endothelial cells",
  "ly Endothelial cells", "Macrophages", "Macrophages M2",
  "CD8+ T-cells", "StromaScore", "ImmuneScore"
)
celulas <- celulas[celulas %in% rownames(xc)]

mat <- xc[celulas, , drop = FALSE]
ann <- data.frame(Group = pheno$group, row.names = pheno$sample)
ann <- ann[colnames(mat), , drop = FALSE]

pheatmap(
  mat,
  scale = "row",
  annotation_col = ann,
  annotation_colors = list(Group = c(CTRL = "#1B9E77",
                                     DCM  = "#D95F02",
                                     ICM  = "#7570B3")),
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize = 10,
  main = "xCell: escores de tipos celulares por amostra",
  filename = file.path(FIGURES, "xCell_heatmap.png"),
  width = 9.5, height = 6.5
)

# Violinos dos tipos selecionados, com ICM como contexto
long <- as.data.frame(t(mat)) |>
  tibble::rownames_to_column("sample") |>
  pivot_longer(-sample, names_to = "cell_type", values_to = "score") |>
  left_join(pheno, by = "sample")

p_violino <- ggplot(long, aes(group, score, fill = group)) +
  geom_violin(width = 0.85, alpha = 0.35, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.14, fill = "white", alpha = 0.9,
               outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.09, size = 0.7, alpha = 0.55, color = "grey20") +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(CTRL = "#1B9E77",
                               DCM  = "#D95F02",
                               ICM  = "#7570B3")) +
  labs(
    title = "xCell: distribuição dos escores por grupo",
    subtitle = paste0("Violino = distribuição dos escores; caixa interna = mediana e quartis. ",
                      "ICM mostrada como contexto; os testes comparam DCM vs CTRL."),
    x = NULL, y = "Escore xCell (enriquecimento relativo)",
    fill = "Grupo"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(hjust = 0, size = 9))

save_plot(p_violino, "xCell_violinos.png", 10, 7)

message("SVA: output/tables/SVA_comparacao.csv e SVA_logFC_antes_depois.png")
message("xCell: output/tables/xCell_teste.csv e xCell_heatmap.png")
step("FIM DO SCRIPT 06 — batch e deconvolução concluídos.")
