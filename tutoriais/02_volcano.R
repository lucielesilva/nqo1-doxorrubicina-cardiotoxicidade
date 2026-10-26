#!/usr/bin/env Rscript

# ==============================================================================
# 02 — VOLCANO PLOT (expressão diferencial DCM vs CTRL)
# ==============================================================================
#
# ESTE É O "FIM DA AULA ATÉ O VOLCANO".
#   Aqui fazemos:
#     1) carregar os dados já preparados;
#     2) comparar DCM vs CTRL com o modelo limma, ajustando por sexo e idade;
#     3) gerar a tabela de expressão diferencial;
#     4) desenhar o volcano plot com o NQO1 destacado.
#
#   Os scripts seguintes (04_rede_PPI, 05_extras_GSEA_GSVA) NÃO são necessários
#   para entender o volcano. Eles continuam a história depois dele.
#
# PERGUNTA BIOLÓGICA:
#   A doxorrubicina (DOX) causa cardiotoxicidade. A enzima NQO1 é citoprotetora
#   e detoxifica quinonas. A hipótese é que, na cardiomiopatia dilatada (DCM),
#   o NQO1 está DOWNREGULADO — ou seja, o coração doente tem MENOS defesa
#   antioxidante e fica mais vulnerável à DOX.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# 1) CARREGAR A BASE
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
  library(limma)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)   # rótulos que "empurram" os textos para não se sobrepor
})

step("Volcano da expressão diferencial (DCM vs CTRL)")

# ------------------------------------------------------------------------------
# 2) DADOS
# ------------------------------------------------------------------------------
dados <- preparar_dados()
expr  <- dados$expr
pheno <- dados$pheno

# ------------------------------------------------------------------------------
# 3) ESCOLHER AS AMOSTRAS DA COMPARAÇÃO
# ------------------------------------------------------------------------------
# Por que excluir a ICM?
#   ICM é outra doença (isquêmica), com biologia diferente da DCM. Misturar as
#   duas poderia criar diferenças que não são da DCM. Para a pergunta do
#   projeto, comparamos apenas coração doente dilatado (DCM) vs coração normal
#   (CTRL). A ICM fica de fora deste teste (mas aparece no UMAP como contexto).
#
# Por que "arrange(match(...))"?
#   Garantimos que a ordem das amostras na tabela de fenótipos seja a MESMA da
#   matriz de expressão. Isso evita o clássico erro de comparar a amostra errada.
pheno_de <- pheno |>
  filter(group %in% c("CTRL", "DCM")) |>
  mutate(group = factor(group, levels = c("CTRL", "DCM"))) |>
  arrange(match(sample, colnames(expr)))

message("Amostras no teste: CTRL=", sum(pheno_de$group == "CTRL"),
        " | DCM=", sum(pheno_de$group == "DCM"))

# Confere se sexo e idade estão completos. O modelo precisa dessas covariáveis.
if (anyNA(pheno_de$sex) || anyNA(pheno_de$age)) {
  stop("Há valores ausentes em sexo ou idade. O modelo não pode ser ajustado.",
       call. = FALSE)
}

# ------------------------------------------------------------------------------
# 4) MATRIZ DE EXPRESSÃO SÓ COM AS AMOSTRAS DO TESTE
# ------------------------------------------------------------------------------
m <- expr[, pheno_de$sample, drop = FALSE]

# Filtro: mantém genes com RPKM > 1 em pelo menos 25% das amostras.
# Genes quase sempre zero não têm informação confiável e só aumentam o número
# de testes (e o problema de múltiplas comparações).
m <- m[rowSums(m > 1) >= 0.25 * ncol(m), , drop = FALSE]

# Transformação log2. Como já explicado: +1 vira dobro, -1 vira metade.
expr_de <- log2(m + 1)
message("Genes no modelo: ", nrow(expr_de))

# ------------------------------------------------------------------------------
# 5) MODELO LIMMA
# ------------------------------------------------------------------------------
# O QUE O limma FAZ?
#   Para CADA gene, ele ajusta um modelo linear:
#
#       expressão ~ grupo (DCM vs CTRL) + sexo + idade
#
#   Traduzindo: "a expressão desse gene difere entre DCM e CTRL, depois de
#   descontar o efeito de sexo e idade?".
#
# POR QUE AJUSTAR POR SEXO E IDADE?
#   Homens/mulheres e pessoas mais velhas/mais novas podem ter expressão
#   diferente por motivos que não são a doença. O ajuste evita confundir
#   "efeito da DCM" com "efeito da idade" ou "efeito do sexo".
#
# O QUE É O DESIGN?
#   É a matriz que diz ao R quem é quem. "~ 0 + group + sex + age" cria uma
#   coluna para CTRL e outra para DCM, sem intercepto, e colunas para sexo e
#   idade. Assim o contraste DCM - CTRL fica direto.
design <- model.matrix(~ 0 + group + sex + age, data = pheno_de)
colnames(design) <- sub("^group", "", colnames(design))

# CONTRATE: é a comparação que queremos. "DCM-CTRL" = DCM menos CTRL.
#   log2FC > 0  -> gene mais expresso na DCM.
#   log2FC < 0  -> gene menos expresso na DCM.
contraste <- makeContrasts("DCM-CTRL", levels = design)

# lmFit ajusta o modelo gene a gene. contrasts.fit aplica o contraste.
# eBayes usa "empirical Bayes": ele "empresta" informação entre genes para
# estabilizar as variâncias, o que ajuda muito quando o n é pequeno.
# trend = TRUE permite que a variância dependa do nível médio de expressão.
fit <- eBayes(contrasts.fit(lmFit(expr_de, design), contraste), trend = TRUE)

# topTable devolve a tabela final. number = Inf = todos os genes.
# sort.by = "none" mantém a ordem original (não ordena por p-valor).
tt <- topTable(fit, number = Inf, sort.by = "none")

# ------------------------------------------------------------------------------
# 6) TABELA DE EXPRESSÃO DIFERENCIAL
# ------------------------------------------------------------------------------
# Colunas:
#   gene            = símbolo do gene
#   log2FoldChange  = tamanho e direção do efeito (DCM vs CTRL)
#   AveExpr         = expressão média do gene
#   pvalue          = p-valor bruto
#   padj            = p-valor ajustado por múltiplos testes (FDR)
#   neg_log10_padj  = -log10(padj), usado no eixo y do volcano
#   direction       = "up" se subiu na DCM, "down" se desceu
#   significant     = TRUE se padj < 0.05
#
# POR QUE padj E NÃO pvalue?
#   Estamos testando ~20 mil genes ao mesmo tempo. Se usarmos p < 0.05 em cada
#   um, muitos genes "significativos" apareceriam só por azar (falso positivo).
#   O padj (Benjamini-Hochberg/FDR) corrige isso: controla a proporção de
#   falsos positivos entre os genes chamados de significativos.
deg <- data.frame(
  gene           = rownames(tt),
  log2FoldChange = tt$logFC,
  AveExpr        = tt$AveExpr,
  pvalue         = tt$P.Value,
  padj           = tt$adj.P.Val,
  stringsAsFactors = FALSE
) |>
  mutate(
    # pmax(padj, 1e-300) evita log10(0) = infinito.
    neg_log10_padj = -log10(pmax(padj, 1e-300)),
    direction      = ifelse(log2FoldChange > 0, "up", "down"),
    significant    = padj < 0.05
  )

write_csv(deg, file.path(TABLES, "expressao_diferencial_completa.csv"))

# ------------------------------------------------------------------------------
# 7) O GENE-ALVO: NQO1
# ------------------------------------------------------------------------------
# Esta é a "estrela" do projeto. Queremos saber se ele está downregulado na DCM.
nqo1 <- filter(deg, gene == "NQO1")
write_csv(nqo1, file.path(TABLES, "expressao_nqo1_diferencial.csv"))

if (nrow(nqo1) == 1) {
  message(sprintf("NQO1: log2FC = %.3f | padj = %.2e | %s",
                  nqo1$log2FoldChange, nqo1$padj, nqo1$direction))
  nqo1_caption <- sprintf(
    "NQO1 destacado: log2FC = %.2f | padj = %.1e (%s na DCM).",
    nqo1$log2FoldChange, nqo1$padj,
    ifelse(nqo1$log2FoldChange < 0, "downregulado", "upregulado"))
} else {
  warning("NQO1 não foi encontrado na tabela de expressão diferencial.")
  nqo1_caption <- "NQO1 não encontrado na matriz de expressão."
}

# ------------------------------------------------------------------------------
# 8) VOLCANO PLOT
# ------------------------------------------------------------------------------
# COMO LER O VOLCANO?
#   * Cada ponto é UM gene.
#   * Eixo x = log2 fold change: tamanho e direção do efeito.
#       direita  = mais expresso na DCM
#       esquerda = menos expresso na DCM
#   * Eixo y = -log10(padj): confiança estatística.
#       mais alto = padj menor = mais significativo.
#   * As linhas tracejadas são os limiares:
#       |log2FC| >= 0.585 (~1,5 vezes) e padj < 0.05.
#   * Genes no canto superior direito = up na DCM e significativos.
#   * Genes no canto superior esquerdo = down na DCM e significativos.
#   * Genes no meio/baixo = pouco efeito ou pouca confiança.
#
# ANALOGIA:
#   O volcano é um "mapa de importância". O eixo horizontal diz "quanto mudou";
#   o vertical diz "quão certeza temos de que mudou". O ideal é olhar os dois
#   juntos: um gene pode mudar muito, mas com pouca confiança, ou mudar pouco,
#   mas de forma muito consistente.
log2fc_cut <- 0.585   # |log2FC| >= 0.585 ~ 1,5x
padj_cut   <- 0.05

volcano <- deg |>
  mutate(regulation = case_when(
    padj < padj_cut & log2FoldChange >=  log2fc_cut ~ "Up na DCM",
    padj < padj_cut & log2FoldChange <= -log2fc_cut ~ "Down na DCM",
    TRUE ~ "Não significativo"
  )) |>
  # factor com níveis definidos para a legenda sair sempre na mesma ordem.
  mutate(regulation = factor(regulation,
                             levels = c("Up na DCM", "Down na DCM",
                                        "Não significativo")))

n_up   <- sum(volcano$regulation == "Up na DCM")
n_down <- sum(volcano$regulation == "Down na DCM")

# Rótulos: top 15 significativos (por padj) + o NQO1, sem repetir.
top15 <- volcano |>
  filter(regulation != "Não significativo") |>
  arrange(padj, desc(abs(log2FoldChange))) |>
  slice_head(n = 15)

labels <- bind_rows(top15, filter(volcano, gene == "NQO1")) |>
  distinct(gene, .keep_all = TRUE)

p_volcano <- ggplot(volcano, aes(log2FoldChange, neg_log10_padj)) +
  # Todos os genes. alpha deixa os pontos mais transparentes para não virar
  # uma "mancha" só; size controla o tamanho.
  geom_point(aes(color = regulation), alpha = 0.6, size = 1.4) +
  # Limiares de efeito (vertical) e de significância (horizontal).
  geom_vline(xintercept = c(-log2fc_cut, log2fc_cut),
             linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(padj_cut),
             linetype = "dashed", color = "grey60") +
  # Círculo extra em volta do NQO1 para ele não passar despercebido.
  geom_point(data = filter(volcano, gene == "NQO1"),
             shape = 21, size = 3.6, color = "black",
             fill = NA, stroke = 0.8) +
  # Rótulos dos top 15 significativos. O ggrepel afasta os textos para evitar
  # que um cubra o outro.
  geom_text_repel(
    data = filter(labels, gene != "NQO1"),
    aes(label = gene),
    size = 2.8, color = "grey20",
    max.overlaps = Inf,
    box.padding = 0.35, point.padding = 0.2,
    min.segment.length = 0, segment.color = "grey50",
    seed = 42
  ) +
  # O NQO1 ganha um rótulo PRÓPRIO, em negrito e com fundo branco. Assim ele
  # aparece mesmo quando há muitos genes significativos ao redor.
  geom_label_repel(
    data = filter(labels, gene == "NQO1"),
    aes(label = gene),
    size = 3.2, fontface = "bold", color = "black", fill = "white",
    label.size = 0.25,
    box.padding = 0.6, point.padding = 0.4,
    min.segment.length = 0, segment.color = "black",
    max.overlaps = Inf, seed = 42
  ) +
  scale_color_manual(values = c("Up na DCM"         = "#1B9E77",
                                "Down na DCM"        = "#D95F02",
                                "Não significativo"  = "grey75")) +
  labs(
    title    = "Volcano da expressão diferencial (DCM vs CTRL)",
    subtitle = paste0("GSE116250 | limma ajustado por sexo e idade | ",
                      "|log2FC| >= 0.585 e padj < 0.05 | ",
                      "Up: ", n_up, " | Down: ", n_down),
    x = "log2 fold change (DCM vs CTRL)",
    y = "-log10(padj)",
    color = "Regulação (DCM vs CTRL)",
    caption = nqo1_caption
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0))

save_plot(p_volcano, "volcano_plot.png", 9, 6.5)

# ------------------------------------------------------------------------------
# 9) SALVAR O RESULTADO DO limma PARA OS PRÓXIMOS SCRIPTS
# ------------------------------------------------------------------------------
# O 03 (PPI) e o 04 (GSEA/GSVA) reaproveitam esta tabela e a matriz log2.
# Assim não precisamos rodar o limma de novo em cada script.
resultado_limma <- list(
  deg       = deg,        # tabela de expressão diferencial completa
  pheno_de  = pheno_de,   # amostras usadas no teste
  expr_de   = expr_de,    # matriz log2 usada no modelo
  design    = design,     # matriz do modelo
  gse       = GSE_ID
)
saveRDS(resultado_limma,
        file.path(PROCESSED, "resultado_limma_DCM_vs_CTRL.rds"))

step("Tabela completa: output/tables/expressao_diferencial_completa.csv")
step("Tabela do NQO1: output/tables/expressao_nqo1_diferencial.csv")
step("Volcano: output/figures/volcano_plot.png")
step("FIM DO SCRIPT 02 — A AULA PODE PARAR AQUI (até o volcano).")
message("Para continuar a história: 03_gene_sets.R -> 04_rede_PPI.R (e, se quiser, 05_extras_GSEA_GSVA.R).")
