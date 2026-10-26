#!/usr/bin/env Rscript

# ==============================================================================
# executar_tutoriais.R
# ------------------------------------------------------------------------------
# Executa, na ordem correta, os scripts didáticos da pasta tutoriais/.
#
# Uso (a partir da raiz do projeto):
#   Rscript executar_tutoriais.R
#
# A ordem de execução é a MESMA da numeração dos arquivos:
#   00 (base) → 01 → 02 → 03 → 04 → 05 → 06
#
# O script 00 é carregado uma única vez aqui, para que os dados do GEO e a via
# AMPK sejam baixados uma só vez e reaproveitados pelos demais scripts.
#
# O script resolve a raiz a partir do próprio arquivo, de modo que funciona
# mesmo quando chamado de outra pasta.
# ==============================================================================

args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- if (length(args)) {
  dirname(normalizePath(sub("^--file=", "", args[1]),
                        winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

setwd(root)
options(stringsAsFactors = FALSE)
set.seed(42)

message("Raiz do projeto: ", root)

# Carrega configuração, caminhos e funções auxiliares.
source(file.path("tutoriais", "00_configuracao_e_dados.R"), encoding = "UTF-8")

# Prepara uma única vez os dados e a via AMPK.
dados <- preparar_dados()
ampk  <- baixar_genes_ampk()

# Ordem de execução = numeração dos arquivos.
scripts <- c(
  "01_UMAP.R",
  "02_volcano.R",
  "03_gene_sets.R",
  "04_rede_PPI.R",
  "05_extras_GSEA_GSVA.R",
  "06_batch_deconvolution.R"
)

faltando <- scripts[!file.exists(file.path("tutoriais", scripts))]
if (length(faltando) > 0) {
  stop("Scripts não encontrados em tutoriais/: ",
       paste(faltando, collapse = ", "), call. = FALSE)
}

for (s in scripts) {
  message("\n", strrep("=", 70))
  message("Executando: ", s)
  message(strrep("=", 70))
  source(file.path("tutoriais", s), encoding = "UTF-8")
}

message("\nPipeline concluído. Resultados em output/ e data/processed/.")
