#!/usr/bin/env Rscript

# ==============================================================================
# 04 — REDE PPI DA VIA AMPK (interações entre proteínas)
# ==============================================================================
#
# ONDE ESTAMOS NA HISTÓRIA?
#   O volcano olhou gene por gene. Agora queremos olhar o SISTEMA: como as
#   proteínas codificadas pelos genes da via AMPK se conectam umas às outras.
#
# O QUE É UMA REDE PPI?
#   PPI = Protein-Protein Interaction (interação proteína-proteína).
#   É um grafo em que:
#     * cada NÓ é uma proteína (representada pelo símbolo do gene);
#     * cada ARESTA (linha) é uma interação física ou funcional entre elas.
#
# ANALOGIA PARA LEIGOS:
#   Pense em uma rede social. Algumas proteínas são "influenciadoras": falam
#   com muitas outras. Essas são os HUBS. Se um hub falha, a rede inteira sente.
#   Na doença, um hub pode ser um ponto frágil ou um alvo terapêutico.
#
# DE ONDE VÊM AS INTERAÇÕES?
#   Do STRING (https://string-db.org), uma base pública que integra evidências
#   de interação. Usamos score >= 700 (escala 0–1000) = alta confiança.
#
# O QUE ESTE SCRIPT FAZ?
#   1) Carrega a tabela do limma gerada pelo 02_volcano.R.
#   2) Baixa a lista de genes da via AMPK do KEGG.
#   3) Recorta a via dentro da tabela de expressão diferencial.
#   4) Baixa as interações do STRING e monta a rede.
#   5) Desenha a rede: cor = log2FC; tamanho = grau (nº de conexões).
#
# IMPORTANTE: rode o 02_volcano.R ANTES deste script.
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
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)   # rótulos que se afastam para não se sobrepor
  library(igraph)
})

step("Rede PPI da via AMPK")

# ------------------------------------------------------------------------------
# 2) CARREGAR O RESULTADO DO limma (produzido pelo script 02)
# ------------------------------------------------------------------------------
rds_limma <- file.path(PROCESSED, "resultado_limma_DCM_vs_CTRL.rds")
if (!file.exists(rds_limma)) {
  stop("Não encontrei o resultado do limma.\n",
       "Rode primeiro o script 02_volcano.R (ele cria o arquivo\n",
       "data/processed/resultado_limma_DCM_vs_CTRL.rds).", call. = FALSE)
}

resultado <- readRDS(rds_limma)
deg <- resultado$deg   # tabela com gene, log2FoldChange, padj etc.

# ------------------------------------------------------------------------------
# 3) GENES DA VIA AMPK
# ------------------------------------------------------------------------------
# baixar_genes_ampk() usa cache: se o CSV já existe, não baixa de novo.
ampk <- baixar_genes_ampk()

# Recortamos da tabela completa só os genes que pertencem à via AMPK.
deg_ampk <- deg |>
  filter(gene %in% ampk$symbol)

write_csv(deg_ampk, file.path(TABLES, "expressao_diferencial_via_ampk.csv"))
message("Genes da via AMPK mapeados na expressão: ", nrow(deg_ampk))

# ------------------------------------------------------------------------------
# 4) OBTER AS INTERAÇÕES DO STRING
# ------------------------------------------------------------------------------
# Tentamos primeiro a API pública do STRING: ela é rápida e não precisa
# baixar as bases grandes do pacote. Se a API falhar (por exemplo, sem
# internet), tentamos o pacote STRINGdb como alternativa.
obter_ppi_string <- function(genes, score_min = 700) {

  # ---- Opção 1: API pública do STRING --------------------------------------
  # A API aceita vários identificadores separados por %0d (quebra de linha).
  ppi <- tryCatch({
    ids <- paste(utils::URLencode(genes, reserved = TRUE), collapse = "%0d")
    url <- paste0(
      "https://string-db.org/api/tsv/network?identifiers=", ids,
      "&species=9606&required_score=", score_min,
      "&caller_identity=nqo1_aula"
    )
    tab <- utils::read.delim(url, stringsAsFactors = FALSE)
    if (nrow(tab) == 0) return(NULL)

    data.frame(
      from  = tab$preferredName_A,
      to    = tab$preferredName_B,
      score = as.numeric(tab$score),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("API do STRING falhou: ", conditionMessage(e))
    NULL
  })

  if (!is.null(ppi) && nrow(ppi) > 0) return(ppi)

  # ---- Opção 2: pacote STRINGdb (fallback) --------------------------------
  # O STRINGdb baixa bases grandes na primeira vez; por isso ele é o plano B.
  if (requireNamespace("STRINGdb", quietly = TRUE)) {
    message("Tentando o pacote STRINGdb como alternativa...")
    ppi <- tryCatch({
      db <- STRINGdb::STRINGdb$new(
        version         = "11.5",
        species         = 9606,          # 9606 = Homo sapiens
        score_threshold = score_min,
        input_directory = STRING_CACHE
      )

      # map() traduz símbolo do gene -> ID interno do STRING.
      mapped <- db$map(data.frame(gene = genes), "gene",
                       removeUnmappedRows = TRUE)

      # get_interactions() baixa as arestas entre os IDs mapeados.
      inter  <- db$get_interactions(mapped$STRING_id)

      # lookup: ID STRING -> símbolo legível (ex.: "9606.ENSP..." -> "PRKAA1").
      lookup <- setNames(mapped$gene, mapped$STRING_id)

      data.frame(
        from  = lookup[inter$from],
        to    = lookup[inter$to],
        score = inter$combined_score,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message("STRINGdb falhou: ", conditionMessage(e))
      NULL
    })

    if (!is.null(ppi) && nrow(ppi) > 0) return(ppi)
  } else {
    message("Pacote STRINGdb não instalado; usei apenas a API pública.")
  }

  NULL
}

ppi <- obter_ppi_string(deg_ampk$gene, score_min = 700)

if (is.null(ppi) || nrow(ppi) == 0) {
  stop("Não foi possível obter a rede PPI do STRING.\n",
       "Verifique a conexão com a internet ou instale o pacote STRINGdb:\n",
       "  BiocManager::install('STRINGdb')", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 5) LIMPAR AS ARESTAS
# ------------------------------------------------------------------------------
# * Remove arestas com gene não mapeado (NA) ou auto-interação (from == to).
# * Remove duplicatas: A-B e B-A são a mesma interação em rede não direcionada.
#   Para isso, criamos uma "chave" com os dois nomes em ordem alfabética.
ppi <- ppi |>
  filter(!is.na(from), !is.na(to), from != to) |>
  filter(from %in% deg_ampk$gene, to %in% deg_ampk$gene) |>
  mutate(par = paste(pmin(from, to), pmax(from, to))) |>
  distinct(par, .keep_all = TRUE) |>
  select(-par)

message("Interações de alta confiança (score >= 700): ", nrow(ppi))

# A API do STRING devolve o score entre 0 e 1; o pacote STRINGdb devolve entre
# 0 e 1000. Padronizamos para 0–1000 para a tabela e para a espessura das linhas.
if (nrow(ppi) > 0 && max(ppi$score, na.rm = TRUE) <= 1) {
  ppi$score <- ppi$score * 1000
}

# ------------------------------------------------------------------------------
# 6) MONTAR O GRAFO
# ------------------------------------------------------------------------------
# graph_from_data_frame() cria o grafo a partir da tabela de arestas.
# directed = FALSE porque interação proteína-proteína não tem direção.
g <- graph_from_data_frame(ppi, directed = FALSE)

# simplify() remove arestas repetidas e laços, caso ainda existam.
# Se duas evidências ligarem o mesmo par, ficamos com o maior score.
g <- simplify(g, edge.attr.comb = list(score = "max", "ignore"))

# Layout "force-directed" (Fruchterman-Reingold): proteínas que interagem se
# atraem; as que não interagem se afastam. É como um sistema de molas.
# niter = 3000 e grid = "nogrid" deixam o desenho mais espalhado e menos
# "embolado" do que o padrão. set.seed(42) no 00 torna o layout reprodutível.
layout <- layout_with_fr(g, weights = E(g)$score, niter = 3000, grid = "nogrid")

# ------------------------------------------------------------------------------
# 7) TABELA DE NÓS (para desenhar)
# ------------------------------------------------------------------------------
# Para cada proteína:
#   x, y           = posição no layout
#   log2FoldChange = efeito na DCM (vem do limma)
#   degree         = grau = número de conexões (quanto maior, mais "hub")
nodes <- data.frame(
  gene = V(g)$name,
  x    = layout[, 1],
  y    = layout[, 2],
  stringsAsFactors = FALSE
) |>
  left_join(select(deg_ampk, gene, log2FoldChange), by = "gene") |>
  mutate(
    degree = as.numeric(degree(g, v = V(g)$name)),
    # Se algum gene não tiver log2FC, usamos 0 (branco no gradiente).
    log2FC_plot = ifelse(is.na(log2FoldChange), 0, log2FoldChange)
  )

# Limite simétrico para o gradiente de cor: o branco fica no zero e as cores
# ficam comparáveis entre genes up e down.
lim_fc <- max(abs(nodes$log2FC_plot), na.rm = TRUE)
if (!is.finite(lim_fc) || lim_fc == 0) lim_fc <- 1
nodes <- nodes |>
  mutate(log2FC_plot = pmax(pmin(log2FC_plot, lim_fc), -lim_fc))

# ------------------------------------------------------------------------------
# 8) TABELA DE ARESTAS (para desenhar as linhas)
# ------------------------------------------------------------------------------
# Precisamos das coordenadas x/y de cada ponta da aresta.
edges <- as_data_frame(g, what = "edges") |>
  left_join(select(nodes, gene, x, y), by = c("from" = "gene")) |>
  rename(x_from = x, y_from = y) |>
  left_join(select(nodes, gene, x, y), by = c("to" = "gene")) |>
  rename(x_to = x, y_to = y)

# Salvamos as tabelas para o aluno explorar no Excel.
write_csv(ppi,   file.path(TABLES, "ppi_ampk_string_arestas.csv"))
write_csv(nodes, file.path(TABLES, "ppi_ampk_string_nos.csv"))

# Tabela só com os hubs (as proteínas mais conectadas) para a aula.
hubs <- nodes |>
  arrange(desc(degree)) |>
  select(gene, degree, log2FoldChange) |>
  slice_head(n = 20)
write_csv(hubs, file.path(TABLES, "ppi_ampk_hubs.csv"))

# ------------------------------------------------------------------------------
# 9) GRÁFICO DA REDE
# ------------------------------------------------------------------------------
# COMO LER A REDE?
#   * Cada bolinha = uma proteína.
#   * Cada linha = uma interação de alta confiança.
#   * Cor: verde = gene up na DCM; laranja = down na DCM; branco = sem dado.
#   * Tamanho: grau. Bolinha grande = hub, conversa com muitas proteínas.
#   * Nós muito conectados são candidatos a "pontos frágeis" da via.
#
# DICA DE AULA:
#   O NQO1 NÃO faz parte da lista de genes da via AMPK do KEGG; por isso ele
#   não aparece nesta rede. A rede mostra o "sistema" da via AMPK; o NQO1 é o
#   alvo do projeto e aparece no volcano. Isso rende uma boa discussão: um gene
#   pode ser importante para a doença sem estar nessa via específica.
p_ppi <- ggplot() +
  # Arestas primeiro (ficam atrás das bolinhas). A espessura representa a
  # confiança da interação no STRING: linha mais grossa = score maior.
  geom_segment(data = edges,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
                   linewidth = score),
               color = "grey70", alpha = 0.35) +
  scale_linewidth_continuous(range = c(0.15, 1.1), guide = "none") +
  # Nós: shape 21 permite cor de preenchimento (fill) + borda.
  geom_point(data = nodes,
             aes(x = x, y = y, size = degree, fill = log2FC_plot),
             shape = 21, color = "grey20", stroke = 0.4) +
  # Se o NQO1 estiver na rede, ele ganha um anel preto. Se não estiver, esta
  # camada simplesmente não desenha nada (o NQO1 não faz parte da via AMPK).
  geom_point(data = filter(nodes, gene == "NQO1"),
             shape = 21, size = 5, color = "black", fill = NA, stroke = 1) +
  scale_fill_gradient2(
    low = "#D95F02", mid = "white", high = "#1B9E77", midpoint = 0,
    limits = c(-lim_fc, lim_fc),
    name = "log2FC\n(DCM vs CTRL)"
  ) +
  scale_size_continuous(
    range = c(2, 9),
    name = "Grau\n(nº de interações)"
  ) +
  labs(
    title    = "Rede PPI da via AMPK na DCM",
    subtitle = sprintf("STRING Homo sapiens, score >= 700 | %d nós e %d arestas | rótulos: 20 hubs",
                       vcount(g), ecount(g)),
    x = NULL, y = NULL,
    caption = paste0("Cor = log2FC (verde up, laranja down na DCM); ",
                     "tamanho = grau; espessura da linha = confiança da interação. ",
                     "Nós com maior grau são hubs da via.")
  ) +
  coord_fixed(clip = "off") +
  theme_void(base_size = 11) +
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0))

# Rótulos: os 20 nós mais conectados. O ggrepel afasta os textos e desenha uma
# linquinha até o nó, deixando a rede muito mais legível.
label_nodes <- nodes |>
  arrange(desc(degree)) |>
  slice_head(n = 20)

# Se o NQO1 estiver na rede, garante que ele também seja rotulado.
if ("NQO1" %in% nodes$gene && !("NQO1" %in% label_nodes$gene)) {
  label_nodes <- bind_rows(label_nodes, filter(nodes, gene == "NQO1"))
}

top5_hubs <- head(hubs$gene, 5)
p_ppi <- p_ppi +
  geom_text_repel(
    data = label_nodes,
    aes(x = x, y = y, label = gene,
        fontface = ifelse(gene %in% top5_hubs, "bold", "plain")),
    size = 3, color = "grey15",
    max.overlaps = Inf,
    box.padding = 0.5, point.padding = 0.3,
    min.segment.length = 0, segment.color = "grey55",
    seed = 42
  )

save_plot(p_ppi, "rede_ppi_ampk.png", 12, 9)

message("Rede PPI: ", vcount(g), " nós e ", ecount(g), " arestas")
message("Top 5 hubs: ",
        paste(head(hubs$gene, 5), collapse = ", "))
step("FIM DO SCRIPT 04 — rede salva em output/figures/rede_ppi_ampk.png")
message("Próximo passo opcional: 05_extras_GSEA_GSVA.R")
