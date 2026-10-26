#!/usr/bin/env Rscript

# ==============================================================================
# 03 — GENES DO DATASET x GENES DA VIA AMPK
# ==============================================================================
# Este script deixa explícita a diferença entre:
#   * todos os genes do dataset após o filtro de expressão;
#   * os genes da via AMPK do KEGG (hsa04152);
#   * a interseção usada no limma, na rede PPI e no GSEA.
#
# A ideia é evitar a confusão "o gene está no dataset?" versus
# "o gene pertence à via AMPK?".
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(42)

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
  library(readr)
  library(ggplot2)
  library(ggforce)
})

dados <- preparar_dados()
ampk  <- baixar_genes_ampk()

rds_limma <- file.path(PROCESSED, "resultado_limma_DCM_vs_CTRL.rds")
if (!file.exists(rds_limma)) {
  stop("Rode primeiro o script 02_volcano.R.", call. = FALSE)
}
deg <- readRDS(rds_limma)$deg

dataset_genes <- deg$gene
ampk_genes    <- ampk$symbol
inter         <- intersect(dataset_genes, ampk_genes)
only_dataset  <- setdiff(dataset_genes, ampk_genes)
only_ampk     <- setdiff(ampk_genes, dataset_genes)

resumo <- data.frame(
  conjunto = c(
    "Dataset após filtro de expressão",
    "Via AMPK (KEGG hsa04152)",
    "Interseção (dataset ∩ via AMPK)",
    "Somente no dataset",
    "Somente na via AMPK"
  ),
  n = c(
    length(dataset_genes),
    length(ampk_genes),
    length(inter),
    length(only_dataset),
    length(only_ampk)
  )
)
write_csv(resumo, file.path(TABLES, "genes_dataset_vs_ampk.csv"))

# ------------------------------------------------------------------------------
# Diagrama de Venn simplificado (dois círculos)
# ------------------------------------------------------------------------------
p_venn <- ggplot() +
  geom_circle(aes(x0 = 0, y0 = 0, r = 1.65),
              fill = "#1B9E77", alpha = 0.22,
              color = "#1B9E77", linewidth = 1.4) +
  geom_circle(aes(x0 = 1.55, y0 = 0, r = 1.15),
              fill = "#D95F02", alpha = 0.22,
              color = "#D95F02", linewidth = 1.4) +
  annotate("text", x = -0.95, y = 1.95,
           label = paste0("Dataset após filtro\n", length(dataset_genes), " genes"),
           fontface = "bold", color = "#146b63", size = 5.2) +
  annotate("text", x = 2.65, y = 1.55,
           label = paste0("Via AMPK (KEGG)\n", length(ampk_genes), " genes"),
           fontface = "bold", color = "#b34700", size = 5.2) +
  annotate("text", x = 0.35, y = 0,
           label = length(inter), size = 9, fontface = "bold") +
  annotate("text", x = -1.05, y = -0.55,
           label = length(only_dataset), size = 5.5) +
  annotate("text", x = 2.35, y = -0.55,
           label = length(only_ampk), size = 5.5) +
  coord_fixed(clip = "off", xlim = c(-2.8, 4.0), ylim = c(-2.2, 2.6)) +
  labs(
    title = "Genes do dataset x genes da via AMPK",
    caption = paste0(
      "O dataset após filtro tem ", length(dataset_genes),
      " genes. A via AMPK do KEGG tem ", length(ampk_genes),
      " genes. A interseção tem ", length(inter),
      " genes; ", length(only_ampk),
      " genes da via ficam fora do dataset filtrado."
    )
  ) +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.caption = element_text(hjust = 0, size = 11))

save_plot(p_venn, "Venn_dataset_AMPK.png", 8.5, 6.2)

message("Venn: output/figures/Venn_dataset_AMPK.png")
message("Tabela: output/tables/genes_dataset_vs_ampk.csv")
step("FIM DO SCRIPT 03 — separação dos genes concluída.")
