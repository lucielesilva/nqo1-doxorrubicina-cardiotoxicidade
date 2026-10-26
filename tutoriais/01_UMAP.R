#!/usr/bin/env Rscript

# ==============================================================================
# 01 — UMAP DAS AMOSTRAS (a "fotografia" não supervisionada dos dados)
# ==============================================================================
#
# QUAL É A IDEIA DESTE SCRIPT?
#   Antes de testar qual gene muda entre DCM e CTRL, queremos olhar os dados
#   "de cima". O UMAP (Uniform Manifold Approximation and Projection) projeta
#   as amostras em dois eixos para que possamos desenhar cada amostra como um
#   ponto e enxergar agrupamentos e sobreposições.
#
# POR QUE UMAP E NÃO PCA?
#   O PCA é uma projeção LINEAR: ele resume a variância em componentes
#   ortogonais. O UMAP é uma projeção NÃO LINEAR: ele tenta preservar a
#   vizinhança local entre as amostras. Em transcriptomas com relações não
#   lineares, o UMAP pode revelar agrupamentos que o PCA suaviza. Em troca,
#   os eixos do UMAP NÃO têm interpretação biológica direta e as distâncias
#   globais devem ser lidas com cautela.
#
# O QUE É IMPORTANTE ENTENDER?
#   * O UMAP é NÃO SUPERVISIONADO: ele NÃO sabe quem é DCM, CTRL ou ICM.
#     As cores são adicionadas DEPOIS, só para interpretar.
#   * Se as amostras de DCM se separam das de CTRL, isso é evidência de que
#     existe sinal biológico entre os grupos.
#   * Se a separação for por outro fator (sexo, idade, lote, outlier),
#     descobrimos isso ANTES de rodar os testes estatísticos.
#
# ORDEM DA AULA:
#   Primeiro script de análise do roteiro. Ele NÃO depende do limma.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# 1) CARREGAR A BASE (configuração, caminhos e função de dados)
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
  library(uwot)      # UMAP
  library(ggplot2)
  library(readr)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# 2) OBTER OS DADOS
# ------------------------------------------------------------------------------
dados <- preparar_dados()
expr  <- dados$expr
pheno <- dados$pheno

step("UMAP das amostras")
message("Amostras disponíveis: ",
        "CTRL=", sum(pheno$group == "CTRL"), " | ",
        "DCM=",  sum(pheno$group == "DCM"),  " | ",
        "ICM=",  sum(pheno$group == "ICM"))

# ------------------------------------------------------------------------------
# 3) PREPARAR A MATRIZ PARA O UMAP
# ------------------------------------------------------------------------------
# 3.1) Transformação log2(RPKM + 1)
#   RPKM é uma escala muito assimétrica. O log2 reduz a influência de valores
#   extremos. O "+1" evita log(0).
expr_log <- log2(expr + 1)

# 3.2) Filtrar genes pouco expressos
#   Mantemos genes com RPKM > 1 em pelo menos 25% das amostras. Genes quase
#   sempre zero adicionam ruído e não ajudam a separar as amostras.
expr_log <- expr_log[rowSums(expr > 1) >= 0.25 * ncol(expr), , drop = FALSE]
message("Genes após filtro de expressão: ", nrow(expr_log))

# 3.3) Escolher os genes mais variáveis
#   O UMAP precisa de genes que MUDAM entre as amostras. Usamos os 5.000 genes
#   com maior variância. Um número menor reduz ruído; um número maior preserva
#   mais informação, mas pode diluir o sinal.
n_genes_umap <- 5000
vars <- apply(expr_log, 1, var)
top  <- names(sort(vars, decreasing = TRUE))[seq_len(min(n_genes_umap, length(vars)))]
message("Genes usados no UMAP (mais variáveis): ", length(top))

# ------------------------------------------------------------------------------
# 4) RODAR O UMAP
# ------------------------------------------------------------------------------
# uwot::umap() espera as AMOSTRAS nas LINHAS. Nossa matriz tem genes nas
# linhas, então transpomos com t().
#
# Parâmetros principais:
#   n_neighbors = 15  -> quantos vizinhos cada amostra considera. Valores
#                        maiores preservam mais estrutura global; menores
#                        destacam estrutura local.
#   min_dist = 0.1    -> quão próximos os pontos podem ficar no desenho.
#   metric = "euclidean" -> distância usada no espaço de genes.
#   scale = TRUE      -> coloca cada gene na mesma "régua" (desvio-padrão 1).
#   seed = 42         -> torna o resultado reprodutível.
#   n_threads = 1     -> evita variação entre execuções paralelas.
set.seed(42)
umap_coords <- uwot::umap(
  X            = t(expr_log[top, , drop = FALSE]),
  n_neighbors  = 15,
  min_dist     = 0.1,
  metric       = "euclidean",
  scale        = TRUE,
  n_components = 2,
  seed         = 42,
  n_threads    = 1,
  verbose      = FALSE
)

# ------------------------------------------------------------------------------
# 5) MONTAR A TABELA DE COORDENADAS
# ------------------------------------------------------------------------------
# Cada linha = uma amostra; UMAP1 e UMAP2 = posição dela no gráfico.
# O grupo é adicionado só para colorir os pontos (o UMAP não usou essa info).
coords <- data.frame(
  sample = rownames(umap_coords),
  UMAP1  = umap_coords[, 1],
  UMAP2  = umap_coords[, 2],
  group  = pheno$group[match(rownames(umap_coords), pheno$sample)],
  stringsAsFactors = FALSE
)
coords$group <- factor(coords$group, levels = c("CTRL", "DCM", "ICM"))

write_csv(coords, file.path(TABLES, "UMAP_coordenadas.csv"))

# Resumo descritivo dos centroides por grupo (útil para a aula).
resumo <- coords |>
  group_by(group) |>
  summarise(
    n           = n(),
    UMAP1_medio = mean(UMAP1),
    UMAP2_medio = mean(UMAP2),
    sd_UMAP1    = sd(UMAP1),
    sd_UMAP2    = sd(UMAP2),
    .groups     = "drop"
  )
write_csv(resumo, file.path(TABLES, "UMAP_resumo_grupos.csv"))

# ------------------------------------------------------------------------------
# 6) GRÁFICO
# ------------------------------------------------------------------------------
# COMO LER O GRÁFICO?
#   * Cada ponto é UMA AMOSTRA (um coração).
#   * Pontos próximos têm perfis de expressão parecidos.
#   * Os eixos UMAP1/UMAP2 NÃO têm unidade biológica; eles são apenas o
#     desenho bidimensional. Não interprete a distância entre eixos como
#     "mais ou menos expressão".
#   * A elipse resume onde 68% das amostras daquele grupo estão.
#   * Se as elipses de CTRL e DCM quase não se sobrepõem, há boa separação.
#   * Se tudo fica misturado, o sinal global é fraco (ou há outro fator
#     dominando, como sexo, idade ou lote).
p_umap <- ggplot(coords, aes(UMAP1, UMAP2, color = group, fill = group)) +
  stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.68) +
  geom_point(size = 2.8, alpha = 0.9) +
  scale_color_manual(values = c(CTRL = "#1B9E77",
                                DCM  = "#D95F02",
                                ICM  = "#7570B3")) +
  scale_fill_manual(values  = c(CTRL = "#1B9E77",
                                DCM  = "#D95F02",
                                ICM  = "#7570B3")) +
  labs(
    title    = "UMAP das amostras de GSE116250",
    subtitle = paste0("Top ", length(top),
                      " genes mais variáveis | análise NÃO supervisionada"),
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Grupo", fill = "Grupo",
    caption = paste0("Cada ponto é uma amostra. As elipses marcam 68% da dispersão do grupo. ",
                     "O UMAP não usa os rótulos DCM/CTRL/ICM para calcular os eixos; ",
                     "as cores entram só na hora de desenhar. ",
                     "Os eixos UMAP não têm unidade biológica direta.")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0))

save_plot(p_umap, "UMAP_amostras.png", 8, 6.5)

message("UMAP salvo em output/figures/UMAP_amostras.png")
step("FIM DO SCRIPT 01 — figura salva em output/figures/UMAP_amostras.png")
message("Próximo passo: 02_volcano.R (a aula pode terminar no volcano).")
