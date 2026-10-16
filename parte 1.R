#!/usr/bin/env Rscript

# ==============================================================================
# PARTE 1 — NQO1 e via AMPK na cardiotoxicidade da doxorrubicina
# ==============================================================================
#
# ------------------------------------------------------------------------------
# PERGUNTA CENTRAL DO PROJETO
# ------------------------------------------------------------------------------
# A doxorrubicina (DOX) causa cardiotoxicidade. A enzima NQO1 (NAD(P)H quinona
# oxidorredutase 1) é citoprotetora e detoxifica quinonas. A hipótese é que, na
# cardiomiopatia dilatada (DCM), o NQO1 está DOWNREGULADO — ou seja, o coração
# doente tem MENOS defesa antioxidante, o que o torna mais vulnerável à DOX.
#
# ------------------------------------------------------------------------------
# FLUXO DESTE SCRIPT (a ordem importa!)
# ------------------------------------------------------------------------------
#   PASSO 1. Baixar os dados do GEO (GSE116250 — RNA-seq de ventrículo esquerdo).
#   PASSO 2. PCA das amostras (exploração/controle de qualidade ANTES de testar).
#   PASSO 3. Baixar os genes da via AMPK (KEGG hsa04152).
#   PASSO 4. Expressão diferencial (limma): DCM vs CTRL, ajustado por sexo/idade.
#   PASSO 5. Volcano plot do transcriptoma completo (com o NQO1 rotulado).
#   PASSO 6. GSEA da via AMPK (ranking de TODOS os genes, sem corte de p-valor).
#   PASSO 7. GSVA da via AMPK (um escore por AMOSTRA, DCM vs CTRL).
#   PASSO 8. Rede PPI dos genes da via AMPK (STRING).
#
# POR QUE ESSA ORDEM?
#   * O PCA vem PRIMEIRO porque é não supervisionado: ele não conhece os grupos.
#     Se as amostras já se separam por DCM/CTRL, há sinal biológico real.
#   * O limma produz (i) a tabela de genes diferencialmente expressos e (ii) um
#     RANKING de log2FC de todos os genes — matéria-prima do GSEA.
#   * GSEA e GSVA vêm depois do limma porque dependem dele. O GSEA não usa corte
#     de p-valor; a GSVA resume a via em um escore por amostra.
#   * A rede PPI fecha o pipeline: sai do "gene isolado" para o "sistema".
#
# ------------------------------------------------------------------------------
# GLOSSÁRIO RÁPIDO (para aula)
# ------------------------------------------------------------------------------
# RPKM .... Reads Per Kilobase per Million mapped reads. Expressão já normalizada
#           pelo comprimento do gene e pela profundidade do sequenciamento.
# log2(x+1)  Transformação que reduz a assimetria/variância e deixa o valor
#           interpretável como "dobro" (ex.: -1 = metade).
# PCA ..... Principal Component Analysis (não supervisionado). PC1 = direção de
#           maior variabilidade; cada ponto é uma amostra.
# limma ... Modelo linear por gene + "empirical Bayes" (encolhe variâncias).
#           Produz log2FC (efeito) e padj (FDR).
# design .. Matriz do modelo: qual é o grupo e quais covariáveis (sexo, idade).
# padj .... p-valor ajustado por comparações múltiplas (Benjamini-Hochberg/FDR).
# Volcano . x = log2FC, y = -log10(padj): mostra efeito x significância.
# GSEA .... Testa se os genes da via se concentram nas pontas do ranking (sem
#           corte). ES = enrichment score; NES = ES normalizado.
# GSVA .... Transforma a via em UM escore POR AMOSTRA (compara DCM vs CTRL).
# PPI ..... Rede de interações entre proteínas (usamos STRING, score >= 700).
#
# DEPENDÊNCIAS: GEOquery, Biobase, limma, dplyr, ggplot2, readr, igraph,
#               KEGGREST, fgsea, GSVA. Opcional: STRINGdb (usa STRINGdb se existir)
#
# COMO RODAR (PowerShell, na pasta do script):
#   & 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' 'parte 1.R'
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(42)   # reprodutibilidade (layout da rede, permutações do GSEA)

# -----------------------------------------------------------------------------
# 0. CONFIGURAÇÃO (IDs, pastas, bibliotecas, utilitários)
# -----------------------------------------------------------------------------
GSE_ID  <- "GSE116250"    # dataset de GEO (RNA-seq, ventrículo esquerdo humano)
AMPK_ID <- "hsa04152"     # via KEGG "AMPK signaling pathway" (hsa = Homo sapiens)

# Descobre a pasta do projeto a partir da localização do próprio script.
# Assim o script roda de qualquer diretório de trabalho.
setwd("C:/Users/oorie/OneDrive/Documentos/TRABALHOS/NQO1")

args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
ROOT <- if (length(args)) normalizePath(dirname(sub("^--file=", "", args[1])), winslash = "/") else getwd()
setwd(ROOT)

# Estrutura de pastas (criadas automaticamente se não existirem).
RAW     <- file.path(ROOT, "data", "raw")        # dados baixados (GEO, cache STRING)
TABLES  <- file.path(ROOT, "output", "tables")   # tabelas .csv
FIGURES <- file.path(ROOT, "output", "figures")  # figuras .png
for (d in c(RAW, TABLES, FIGURES)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Pacotes obrigatórias. suppressPackageStartupMessages só esconde os avisos
# de carregamento (ex.: funções do dplyr que "mascaram" as do stats).
suppressPackageStartupMessages({
  library(GEOquery); library(Biobase); library(limma); library(KEGGREST)
  library(fgsea);    library(GSVA);    library(igraph)
  library(dplyr);    library(ggplot2); library(readr)
})

# Salva um gráfico com configuração padrão de publicação (300 dpi, fundo branco).
save_plot <- function(p, file, w = 9, h = 6) {
  ggsave(file.path(FIGURES, file), p, width = w, height = h, dpi = 300, bg = "white")
}
# Marcador de passo no console (deixa o log legível).
step <- function(...) message("== ", ...)

# =============================================================================
# PASSO 1. DADOS DO GEO (matriz RPKM + metadados)
# =============================================================================
# O que é o GEO?  Repositório público de dados de expressão do NCBI. Cada estudo
# tem um ID GSE. Baixamos duas coisas:
#   1) SERIES MATRIX -> metadados das amostras (doença, sexo, idade).
#   2) ARQUIVO SUPLEMENTAR (RPKM) -> matriz gene x amostra.
# Sobre o GSE116250: RNA-seq de ventrículo esquerdo (Illumina HiSeq 2500), em
# RPKM (~57.974 genes). Coorte: 14 controles, 37 DCM e 13 ICM. No contraste
# principal excluímos a ICM para comparar DCM vs coração normal.
# =============================================================================
step("Baixando ", GSE_ID, " do GEO")

eset  <- getGEO(GSE_ID, GSEMatrix = TRUE, getGPL = FALSE, destdir = RAW)[[1]]  # ExpressionSet
meta  <- pData(eset)                                                            # tabela de fenótipos
pheno <- data.frame(
  sample  = as.character(meta$title),          # nome da amostra (NF1, DCM8...)
  sex     = as.character(meta[["Sex:ch1"]]),   # covariável
  age     = suppressWarnings(as.numeric(meta[["age:ch1"]])),  # covariável
  stringsAsFactors = FALSE
)

# Traduz os rótulos de doença para três grupos: CTRL (não falência), DCM e ICM.
pheno$group <- dplyr::case_when(
  meta[["disease:ch1"]] == "non-failing"              ~ "CTRL",
  meta[["disease:ch1"]] == "dilated cardiomyopathy"   ~ "DCM",
  meta[["disease:ch1"]] == "ischemic cardiomyopathy"  ~ "ICM",
  TRUE ~ NA_character_
)

pheno <- pheno[!is.na(pheno$group), ]

# Arquivo suplementar de RPKM (se já existir no cache, NÃO baixa de novo).
rpkm_path <- file.path(RAW, paste0(GSE_ID, "_rpkm.txt.gz"))
if (!file.exists(rpkm_path)) {
  getGEOSuppFiles(GSE_ID, baseDir = RAW, fetch_files = TRUE)
  f <- list.files(file.path(RAW, GSE_ID), pattern = "rpkm.*(txt|gz)$", full.names = TRUE)[1]
  file.copy(f, rpkm_path, overwrite = TRUE)
}

# col_types força Gene/Common_name como texto e o resto como número — evita que
# símbolos como "1-Mar" sejam lidos como data.
rpkm  <- as.data.frame(read_tsv(
  rpkm_path, show_col_types = FALSE,
  col_types = cols(.default = col_double(), Gene = "c", Common_name = "c")))
expr  <- as.matrix(rpkm[, -(1:2)])            # remove as 2 colunas de anotação
rownames(expr) <- make.unique(rpkm$Common_name)  # símbolos duplicados ficam ".1"
storage.mode(expr) <- "double"                # garante tipo numérico

# Alinha metadados e matriz: mantém só amostras presentes nos dois.
keep  <- intersect(pheno$sample, colnames(expr))
expr  <- expr[, keep, drop = FALSE]
pheno <- pheno[match(keep, pheno$sample), ]; rownames(pheno) <- NULL
write_csv(pheno, file.path(TABLES, "metadados_das_amostras.csv"))
message("Matriz: ", nrow(expr), " genes x ", ncol(expr), " amostras | CTRL=",
        sum(pheno$group == "CTRL"), " DCM=", sum(pheno$group == "DCM"),
        " ICM=", sum(pheno$group == "ICM"))

# =============================================================================
# PASSO 2. PCA DAS AMOSTRAS (exploração e controle de qualidade)
# =============================================================================
# POR QUE FAZER ISSO PRIMEIRO?
#   O PCA é NÃO SUPERVISIONADO: ele não usa os rótulos DCM/CTRL. Se as amostras
#   se separam por grupo, é evidência de sinal biológico real. Se se separam por
#   outro padrão (lote, sexo, outlier), descobrimos isso ANTES de testar.
# COMO FUNCIONA?
#   1) log2(RPKM+1);  2) filtra genes pouco expressos;  3) seleciona os ~2.000
#   genes mais variáveis;  4) centraliza/escala e roda prcomp();  5) projeta PC1 x PC2.
# =============================================================================
step("PCA das amostras")

expr_log <- log2(expr + 1)   # estabiliza a variância
expr_log <- expr_log[rowSums(expr > 1) >= 0.25 * ncol(expr), ]  # remove genes pouco expressos

# Variância por gene: genes quase constantes só adicionam ruído e distorcem distâncias.
vars   <- apply(expr_log, 1, var)
top    <- names(sort(vars, decreasing = TRUE))[seq_len(min(5000, length(vars)))]
# prcomp espera amostras nas LINHAS -> transpomos. center/scale padronizam cada gene.
pca    <- prcomp(t(expr_log[top, , drop = FALSE]), center = TRUE, scale. = TRUE)
varexp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)   # % de variância por PC

# Coordenadas das amostras no espaço PC1 x PC2 (salvas para inspeção).
coords <- data.frame(
  sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  group  = pheno$group[match(rownames(pca$x), pheno$sample)])
write_csv(coords, file.path(TABLES, "PCA_coordenadas.csv"))

# Gráfico: cada ponto = 1 amostra; a elipse resume a dispersão (68%) do grupo.
save_plot(
  ggplot(coords, aes(PC1, PC2, color = group, fill = group)) +
    stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.68) +
    geom_point(size = 2.6, alpha = 0.9) +
    scale_color_manual(values = c(CTRL = "#1B9E77", DCM = "#D95F02", ICM = "#7570B3")) +
    scale_fill_manual(values  = c(CTRL = "#1B9E77", DCM = "#D95F02", ICM = "#7570B3")) +
    labs(title = "PCA das amostras de GSE116250",
         subtitle = paste("Top", length(top), "genes mais variáveis (não supervisionado)"),
         x = paste0("PC1 (", varexp[1], "% da variância)"),
         y = paste0("PC2 (", varexp[2], "% da variância)"),
         color = "Grupo", fill = "Grupo",
         caption = "Cada ponto é uma amostra; as elipses marcam 68% da dispersão de cada grupo. O PCA não usa os rótulos.") +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom", plot.caption = element_text(hjust = 0)) ,
  "PCA_amostras.png", 8, 6.5)
message("PC1 = ", varexp[1], "% | PC2 = ", varexp[2], "%")

# =============================================================================
# PASSO 3. GENES DA VIA AMPK (KEGG hsa04152)
# =============================================================================
# O que é o KEGG?  Base de vias. A via "AMPK signaling pathway" (hsa04152) é o
# eixo metabólico ligado a estresse energético, autofagia e sobrevivência
# celular no coração. Testamos se a via COMO UM TODO está perturbada na DCM.
# ATENÇÃO: keggGet devolve o campo GENE "plano": posições ímpares = ID Entrez,
# posições pares = "SYMBOL; descrição". Extraímos o símbolo (texto antes do ";").
# =============================================================================
step("Baixando a via AMPK (", AMPK_ID, ") do KEGG")

ampk_path <- file.path(TABLES, "genes_da_via_ampk_kegg.csv")
kegg <- tryCatch(keggGet(AMPK_ID)[[1]]$GENE, error = function(e) NULL)  # tolera falha de rede

if (!is.null(kegg)) {
  defs <- kegg[seq(2, length(kegg), 2)]   # posições pares = "SYMBOL; descrição"
  ampk <- data.frame(
    entrez_id   = kegg[seq(1, length(kegg), 2)],           # posições ímpares = Entrez
    symbol      = trimws(sub(";.*$", "", defs)),           # texto antes do ";"
    description = trimws(defs),
    stringsAsFactors = FALSE)
  ampk <- ampk[ampk$symbol != "", ]                        # remove símbolos vazios
  ampk <- ampk[!duplicated(ampk$symbol), ]                 # remove duplicados
  write_csv(ampk, ampk_path)
} else {
  ampk <- read_csv(ampk_path, show_col_types = FALSE)      # fallback: CSV local
}
message("Via AMPK: ", nrow(ampk), " genes")

# =============================================================================
# PASSO 4. EXPRESSÃO DIFERENCIAL (limma) — DCM vs CTRL
# =============================================================================
# ESTRATÉGIA: ajustamos o modelo no TRANSCRIPTOMA COMPLETO e só depois recortamos
# a via AMPK. Por quê? (i) a correção de múltiplos testes (FDR) precisa de todos
# os genes; (ii) o GSEA exige o ranking de TODOS os genes.
#
# O QUE O limma FAZ?  Para CADA gene, ajusta:
#        expressão ~ grupo (DCM vs CTRL) + sexo + idade
# O coeficiente do grupo é o log2 fold change (log2FC). O limma usa "empirical
# Bayes" para encolher variâncias em direção a um valor comum (estabiliza com n
# pequeno) e aplica FDR (Benjamini-Hochberg) gerando o padj.
# A ICM é EXCLUÍDA — não queremos misturar etiologias.
# =============================================================================
step("Expressão diferencial (limma)")

pheno_de <- pheno |>
  filter(group %in% c("CTRL", "DCM")) |>
  mutate(group = factor(group, levels = c("CTRL", "DCM"))) |>  # CTRL = referência
  arrange(match(sample, colnames(expr)))                        # casa com a matriz

m <- expr[, pheno_de$sample, drop = FALSE]
m <- m[rowSums(m > 1) >= 0.25 * ncol(m), ]   # mantém genes expressos (>1 RPKM em >=25%)
expr_de <- log2(m + 1)                        # transformação log2

# design ~0 + ... cria uma coluna por grupo (sem intercepto); o contraste DCM-CTRL
# fica direto. Sexo e idade entram como covariáveis (ajuste).
design <- model.matrix(~ 0 + group + sex + age, data = pheno_de)
colnames(design) <- sub("^group", "", colnames(design))

fit <- eBayes(contrasts.fit(lmFit(expr_de, design),
                            makeContrasts("DCM-CTRL", levels = design)), trend = TRUE)
tt  <- topTable(fit, number = Inf, sort.by = "none")   # tabela completa (sem ordenar)

deg <- data.frame(
  gene = rownames(tt), log2FoldChange = tt$logFC, AveExpr = tt$AveExpr,
  pvalue = tt$P.Value, padj = tt$adj.P.Val, stringsAsFactors = FALSE) |>
  mutate(
    neg_log10_padj = -log10(pmax(padj, 1e-300)),   # eixo y do volcano (evita log de 0)
    direction      = ifelse(log2FoldChange > 0, "up", "down"),
    significant    = padj < 0.05)
write_csv(deg, file.path(TABLES, "expressao_diferencial_completa.csv"))

# Gene alvo do projeto: NQO1.
nqo1 <- filter(deg, gene == "NQO1")
write_csv(nqo1, file.path(TABLES, "expressao_nqo1_diferencial.csv"))
message(sprintf("NQO1: log2FC = %.3f | padj = %.2e | %s",
                nqo1$log2FoldChange, nqo1$padj, nqo1$direction))

# Recorte da via AMPK a partir do resultado do modelo completo.
deg_ampk <- filter(deg, gene %in% ampk$symbol)
write_csv(deg_ampk, file.path(TABLES, "expressao_diferencial_via_ampk.csv"))
message("Genes da via AMPK mapeados: ", nrow(deg_ampk))

# =============================================================================
# PASSO 5. VOLCANO (transcriptoma completo), com o NQO1 rotulado
# =============================================================================
# COMO LER O VOLCANO?
#   * x = log2 fold change -> magnitude e direção (direita = mais expresso na DCM).
#   * y = -log10(padj)     -> significância (mais alto = menor padj).
#   * limiares: |log2FC| >= 0.585 (~1,5x) e padj < 0.05.
# Rótulos: top 15 significativos + o NQO1 (gene-alvo do projeto).
# =============================================================================
step("Volcano da expressão diferencial")

log2fc_cut <- 0.585   # |log2FC| >= 0.585 ~ 1.5x
padj_cut   <- 0.05

volcano <- deg |>
  mutate(regulation = case_when(
    padj < padj_cut & log2FoldChange >=  log2fc_cut ~ "Up na DCM",
    padj < padj_cut & log2FoldChange <= -log2fc_cut ~ "Down na DCM",
    TRUE ~ "Não significativo"))
n_up   <- sum(volcano$regulation == "Up na DCM")
n_down <- sum(volcano$regulation == "Down na DCM")

# Top 15 significativos (por padj) + o NQO1, sem repetir genes.
top15  <- volcano |>
  filter(regulation != "Não significativo") |>
  arrange(padj, desc(abs(log2FoldChange))) |>
  slice_head(n = 15)
labels <- bind_rows(top15, filter(volcano, gene == "NQO1")) |>
  distinct(gene, .keep_all = TRUE)

save_plot(
  ggplot(volcano, aes(log2FoldChange, neg_log10_padj)) +
    geom_point(aes(color = regulation), alpha = 0.6, size = 1.4) +
    geom_vline(xintercept = c(-log2fc_cut, log2fc_cut), linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "grey60") +
    geom_text(data = labels, aes(label = gene), size = 2.8,
              check_overlap = TRUE, vjust = -0.7, color = "grey15") +
    scale_color_manual(values = c("Up na DCM" = "#1B9E77",
                                  "Down na DCM" = "#D95F02",
                                  "Não significativo" = "grey75")) +
    labs(title = "Volcano da expressão diferencial (DCM vs CTRL)",
         subtitle = paste0("GSE116250 | limma ajustado por sexo e idade | ",
                           "|log2FC| >= 0.585 e padj < 0.05 | ",
                           "Up: ", n_up, " | Down: ", n_down),
         x = "log2 fold change (DCM vs CTRL)",
         y = "-log10(padj)",
         color = "Regulação (DCM vs CTRL)",
         caption = sprintf("NQO1 destacado: log2FC = %.2f | padj = %.1e (downregulado na DCM).",
                           nqo1$log2FoldChange, nqo1$padj)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.caption = element_text(hjust = 0)),
  "volcano_plot.png", 9, 6.5)

# =============================================================================
# PASSO 6. GSEA DA VIA AMPK (fgsea)
# =============================================================================
# POR QUE GSEA E NÃO SÓ A LISTA DE DEGs?  Um corte de padj < 0.05 é arbitrário.
# O GSEA NÃO usa corte: ordena TODOS os genes por log2FC e pergunta se os genes
# da via se concentram nas pontas do ranking.
# CONCEITOS: ES (enrichment score) > 0 = via no topo (up na DCM); ES < 0 = base
# (down). NES = ES normalizado pelo tamanho da via. pval = significância.
# =============================================================================
step("GSEA da via AMPK")

# Ranking: gene -> log2FC, ordenado do maior para o menor (define topo/base).
ranked    <- deg |> filter(!is.na(log2FoldChange)) |> arrange(desc(log2FoldChange))
stats_vec <- setNames(ranked$log2FoldChange, ranked$gene)

gsea <- fgseaSimple(list(AMPK = ampk$symbol), stats_vec, nperm = 10000)[1, ]
# Guardamos o leadingEdge como texto (CSV não aceita coluna-lista).
gsea_out <- data.frame(
  pathway = "KEGG_AMPK_SIGNALING_PATHWAY_HSA04152",
  pval = gsea$pval, padj = gsea$padj, ES = gsea$ES, NES = gsea$NES, size = gsea$size,
  leadingEdge = paste(gsea$leadingEdge[[1]], collapse = "; "))
write_csv(gsea_out, file.path(TABLES, "gsea_via_ampk.csv"))
message(sprintf("GSEA AMPK: NES = %.3f | pval = %.4f", gsea$NES, gsea$pval))

# Curva de enriquecimento: o "random walk" cujo máximo é o ES.
save_plot(
  plotEnrichment(ampk$symbol, stats_vec) +
    labs(title = "GSEA da via AMPK (hsa04152)",
         subtitle = sprintf("NES = %.2f | pval = %.3f | tamanho da via = %d genes",
                            gsea$NES, gsea$pval, gsea$size),
         x = "Ranking dos genes por log2FC (down na DCM <- -> up na DCM)",
         y = "Enrichment score (ES)",
         caption = "Cada traço vertical é um gene da via. ES > 0 = via concentrada no topo (up); ES < 0 = na base (down).") +
    theme_minimal(base_size = 12) +
    theme(plot.caption = element_text(hjust = 0)),
  "GSEA_AMPK_enrichment_plot.png", 8, 6)

# =============================================================================
# PASSO 7. GSVA DA VIA AMPK (escore por amostra)
# =============================================================================
# O QUE MUDA EM RELAÇÃO AO GSEA?  O GSEA dá UM valor para o grupo; a GSVA dá UM
# valor POR AMOSTRA, transformando a via em um "supergene" contínuo. Assim
# comparamos DCM vs CTRL amostra a amostra (boxplot + Wilcoxon).
# kcdf = "Gaussian" porque os dados já são contínuos em log2 (típico de RNA-seq).
# =============================================================================
step("GSVA da via AMPK")

gsva_scores <- gsva(gsvaParam(as.matrix(expr_de), list(AMPK = ampk$symbol), kcdf = "Gaussian"))
gsva_df <- data.frame(
  sample = colnames(gsva_scores),
  AMPK_score = as.numeric(gsva_scores["AMPK", ]),
  group = pheno_de$group[match(colnames(gsva_scores), pheno_de$sample)])
write_csv(gsva_df, file.path(TABLES, "gsva_via_ampk.csv"))

# Teste não paramétrico de Wilcoxon: compara os escores DCM vs CTRL.
p_wilcox <- wilcox.test(AMPK_score ~ group, data = gsva_df)$p.value
message(sprintf("GSVA Wilcoxon DCM vs CTRL: p = %.4f", p_wilcox))

save_plot(
  ggplot(gsva_df, aes(group, AMPK_score, fill = group)) +
    geom_boxplot(width = 0.55, alpha = 0.75, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 1.4, alpha = 0.7, color = "grey25") +
    scale_fill_manual(values = c(CTRL = "#1B9E77", DCM = "#D95F02")) +
    labs(title = "GSVA da via AMPK (escore por amostra)",
         subtitle = sprintf("Wilcoxon DCM vs CTRL: p = %.3f | n = %d amostras", p_wilcox, nrow(gsva_df)),
         x = NULL, y = "Escore GSVA da via AMPK",
         fill = "Grupo",
         caption = "Cada ponto é uma amostra; a caixa mostra a mediana e os quartis. Escore alto = via mais ativa naquela amostra.") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.caption = element_text(hjust = 0)),
  "GSVA_AMPK_boxplot.png", 7.5, 6)

# =============================================================================
# PASSO 8. REDE PPI DA VIA AMPK (STRING, combined_score >= 700)
# =============================================================================
# O QUE É UMA REDE PPI?  Grafo em que nós são proteínas e arestas são interações
# (físicas/funcionais). As interações vêm do STRING (combined_score >= 700 =
# alta confiança). Fecha o raciocínio: a expressão olha gene a gene; a rede
# mostra o SISTEMA. Proteínas muito conectadas (hubs) são pontos frágeis.
# A cor do nó = log2FC (verde = up na DCM, laranja = down); o tamanho = grau.
# =============================================================================
step("Rede PPI (STRING)")

ppi <- NULL
if (requireNamespace("STRINGdb", quietly = TRUE)) {
  ppi <- tryCatch({
    db     <- STRINGdb::STRINGdb$new(version = "11.5", species = 9606,
                                     score_threshold = 700,
                                     input_directory = file.path(RAW, "string_cache"))
    mapped <- db$map(data.frame(gene = deg_ampk$gene), "gene", removeUnmappedRows = TRUE)
    inter  <- db$get_interactions(mapped$STRING_id)
    lookup <- setNames(mapped$gene, mapped$STRING_id)   # ID STRING -> símbolo
    data.frame(from = lookup[inter$from], to = lookup[inter$to],
               score = inter$combined_score) |>
      filter(!is.na(from), !is.na(to), from != to) |>
      distinct(from, to, .keep_all = TRUE)
  }, error = function(e) { message("STRINGdb falhou: ", conditionMessage(e)); NULL })
}

if (!is.null(ppi) && nrow(ppi) > 0) {
  g      <- graph_from_data_frame(ppi, directed = FALSE)   # grafo não direcionado
  layout <- layout_with_fr(g, weights = E(g)$score)        # força-dirigido (Fruchterman-Reingold)
  nodes  <- data.frame(gene = V(g)$name, x = layout[, 1], y = layout[, 2]) |>
    left_join(select(deg_ampk, gene, log2FoldChange), by = "gene") |>
    mutate(degree = as.numeric(igraph::degree(g)))         # grau = nº de conexões
  edges  <- igraph::as_data_frame(g, what = "edges") |>    # igraph:: evita conflito com dplyr
    left_join(select(nodes, gene, x, y), by = c("from" = "gene")) |>
    rename(x_from = x, y_from = y) |>
    left_join(select(nodes, gene, x, y), by = c("to" = "gene")) |>
    rename(x_to = x, y_to = y)

  write_csv(ppi,   file.path(TABLES, "ppi_ampk_string_arestas.csv"))
  write_csv(nodes, file.path(TABLES, "ppi_ampk_string_nos.csv"))

  save_plot(
    ggplot() +
      geom_segment(data = edges, aes(x_from, y_from, xend = x_to, yend = y_to),
                   color = "grey72", alpha = 0.55, linewidth = 0.3) +
      geom_point(data = nodes, aes(x, y, size = degree, fill = log2FoldChange),
                 shape = 21, color = "grey20", stroke = 0.35) +
      geom_text(data = nodes |> arrange(desc(degree)) |> slice_head(n = 35),
                aes(x, y, label = gene), size = 2.8, check_overlap = TRUE, vjust = -0.8) +
      scale_fill_gradient2(low = "#D95F02", mid = "white", high = "#1B9E77") +
      labs(title = "Rede PPI da via AMPK na DCM",
           subtitle = sprintf("STRING Homo sapiens, score >= 700 | %d nós e %d arestas",
                              vcount(g), ecount(g)),
           x = NULL, y = NULL, fill = "log2FC (DCM vs CTRL)", size = "Grau (nº de interações)",
           caption = "Cor = log2FC (verde up, laranja down na DCM); tamanho = grau. Nós com maior grau são hubs da via.") +
      theme_void(base_size = 11) +
      theme(legend.position = "bottom", plot.caption = element_text(hjust = 0)),
    "rede_ppi_ampk.png", 10, 8)
  message("Rede PPI: ", vcount(g), " nós e ", ecount(g), " arestas")
}

step("Pipeline concluído.")
