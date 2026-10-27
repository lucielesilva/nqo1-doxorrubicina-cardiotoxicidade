#!/usr/bin/env Rscript

# ==============================================================================
# 07 — REDE PPI CENTRADA EM NQO1 + VIAS KEGG E REACTOME (versão em R)
# ==============================================================================
#
# O QUE ESTE SCRIPT FAZ?
#   1) Baixa do STRING a rede de vizinhança direta do NQO1 (Homo sapiens,
#      score combinado >= 700, ou seja, alta confiança).
#   2) Consulta o KEGG para obter as vias que contêm NQO1 e os genes dessas vias.
#   3) Consulta o Reactome para obter as vias que contêm NQO1 e sua hierarquia.
#   4) Cruza os parceiros do STRING com as vias KEGG/Reactome.
#   5) Anota os nós com o log2FoldChange (DCM vs CTRL) do resultado do limma.
#   6) Gera figuras (rede, barras, heatmap, mapas de via) e tabelas CSV.
#
# A interpretação dos resultados está integrada ao RELATORIO_CIENTIFICO.md
# (seção 4). Este script é opcional e não faz parte do pipeline principal em R.
#
# COMO RODAR
#   Rscript tutoriais/07_PPI_NQO1_string_kegg_reactome.R
#
# Para forçar a atualização dos dados baixados:
#   NQO1_REFRESH=1 Rscript tutoriais/07_PPI_NQO1_string_kegg_reactome.R
#
# SAÍDAS PRINCIPAIS
#   output/figures/nqo1_rede_ppi_string.png
#   output/figures/nqo1_kegg_pathways.png
#   output/figures/nqo1_reactome_pathways.png
#   output/figures/nqo1_heatmap_vias_parceiros.png
#   output/figures/nqo1_kegg_diagram_hsa05208.png
#   output/figures/nqo1_reactome_diagram_nfe2l2.png
#   output/tables/nqo1_string_ppi_arestas.csv
#   output/tables/nqo1_string_ppi_nos.csv
#   output/tables/nqo1_kegg_pathways.csv
#   output/tables/nqo1_reactome_pathways.csv
#   output/tables/nqo1_reactome_hierarquia.csv
#   output/tables/nqo1_resumo_metricas.csv
#
# OBSERVAÇÃO
#   KEGG e Reactome são bases de vias, não de interações físicas. A rede PPI
#   em si vem do STRING; KEGG/Reactome são usados para anotar e contextualizar
#   os parceiros.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(42)

# ------------------------------------------------------------------------------
# 0) CAMINHOS E CONFIGURAÇÃO
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

detectar_raiz <- function(pasta_script) {
  candidatos <- unique(c(
    pasta_script,
    dirname(pasta_script),
    dirname(dirname(pasta_script))
  ))
  for (d in candidatos) {
    if (!dir.exists(d)) next
    if (file.exists(file.path(d, "DESCRIPTION"))) return(d)
    if (dir.exists(file.path(d, "data")) && dir.exists(file.path(d, "output"))) return(d)
  }
  pasta_script
}

SCRIPT_DIR <- pasta_do_script()
ROOT <- detectar_raiz(SCRIPT_DIR)
setwd(ROOT)

OUT   <- file.path(ROOT, "output")
FIG   <- file.path(OUT, "figures")
TAB   <- file.path(OUT, "tables")
TMP   <- file.path(OUT, "temp")
CACHE <- file.path(TMP, "nqo1_ppi_multidb")
for (d in c(OUT, FIG, TAB, TMP, CACHE)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

GENE          <- "NQO1"
UNIPROT       <- "P15559"
SPECIES       <- 9606
STRING_SCORE  <- 700
STRING_LIMIT  <- 500
REFRESH       <- identical(Sys.getenv("NQO1_REFRESH", "0"), "1")
ALIAS         <- c("SOD2-2" = "SOD2")

log_msg <- function(...) {
  message("[NQO1] ", paste0(..., collapse = ""))
}

alias_gene <- function(x) {
  ifelse(x %in% names(ALIAS), unname(ALIAS[x]), x)
}

# ------------------------------------------------------------------------------
# 1) PACOTES
# ------------------------------------------------------------------------------
pacotes <- c("httr", "jsonlite", "dplyr", "readr", "stringr",
             "igraph", "ggplot2", "magick")
faltando <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltando) > 0) {
  stop("Instale os pacotes: ", paste(faltando, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(readr)
  library(stringr)
  library(igraph)
  library(ggplot2)
  library(magick)
})

# ------------------------------------------------------------------------------
# 2) DOWNLOAD COM CACHE
# ------------------------------------------------------------------------------
get_retry <- function(url, query = list(), type = c("text", "json", "raw")) {
  type <- match.arg(type)
  resposta <- NULL
  ultimo_erro <- NULL
  for (tentativa in seq_len(3)) {
    resposta <- tryCatch({
      if (length(query) > 0) {
        r <- httr::GET(
          url,
          query = query,
          httr::user_agent("nqo1-ppi-analysis/1.0 (educational)"),
          httr::timeout(60)
        )
      } else {
        r <- httr::GET(
          url,
          httr::user_agent("nqo1-ppi-analysis/1.0 (educational)"),
          httr::timeout(60)
        )
      }
      httr::stop_for_status(r)
      r
    }, error = function(e) {
      ultimo_erro <<- e
      NULL
    })
    if (!is.null(resposta)) break
    if (tentativa < 3) Sys.sleep(1.5 * tentativa)
  }
  if (is.null(resposta)) {
    stop("Falha ao acessar ", url, ": ", conditionMessage(ultimo_erro), call. = FALSE)
  }

  switch(
    type,
    text = httr::content(resposta, as = "text", encoding = "UTF-8"),
    json = jsonlite::fromJSON(
      httr::content(resposta, as = "text", encoding = "UTF-8"),
      simplifyVector = FALSE
    ),
    raw = httr::content(resposta, as = "raw")
  )
}

cached_text <- function(nome, url, query = list()) {
  caminho <- file.path(CACHE, nome)
  normalizar <- function(x) gsub("\r\n?", "\n", x)
  if (file.exists(caminho) && !REFRESH) {
    return(normalizar(readr::read_file(caminho)))
  }
  texto <- normalizar(get_retry(url, query = query, type = "text"))
  readr::write_file(texto, caminho)
  texto
}

cached_json <- function(nome, url, query = list()) {
  caminho <- file.path(CACHE, nome)
  if (file.exists(caminho) && !REFRESH) {
    return(jsonlite::fromJSON(caminho, simplifyVector = FALSE))
  }
  dados <- get_retry(url, query = query, type = "json")
  jsonlite::write_json(dados, caminho, auto_unbox = TRUE,
                       pretty = TRUE, null = "null")
  dados
}

cached_raw <- function(nome, url, query = list()) {
  caminho <- file.path(CACHE, nome)
  if (file.exists(caminho) && !REFRESH) {
    return(readBin(caminho, what = "raw", n = as.integer(file.info(caminho)$size)))
  }
  bytes <- get_retry(url, query = query, type = "raw")
  writeBin(bytes, caminho)
  bytes
}

# ------------------------------------------------------------------------------
# 3) STRING — REDE CENTRADA EM NQO1
# ------------------------------------------------------------------------------
baixar_rede_string <- function() {
  log_msg("Baixando rede do STRING...")
  texto <- cached_text(
    "string_nqo1_network.tsv",
    "https://string-db.org/api/tsv/network",
    query = list(
      identifiers = GENE,
      species = SPECIES,
      required_score = STRING_SCORE,
      caller_identity = "nqo1_ppi_multidb",
      limit = STRING_LIMIT
    )
  )

  df <- readr::read_tsv(I(texto), show_col_types = FALSE, progress = FALSE)
  if (nrow(df) == 0) stop("O STRING não retornou interações para NQO1.", call. = FALSE)

  df <- df |>
    dplyr::rename(
      gene_a = preferredName_A,
      gene_b = preferredName_B,
      string_id_a = stringId_A,
      string_id_b = stringId_B
    ) |>
    dplyr::mutate(
      par = apply(cbind(gene_a, gene_b), 1, function(x) paste(sort(x), collapse = "|"))
    ) |>
    dplyr::arrange(dplyr::desc(score)) |>
    dplyr::distinct(par, .keep_all = TRUE) |>
    dplyr::select(-par) |>
    dplyr::mutate(
      is_nqo1_edge = gene_a == GENE | gene_b == GENE,
      score_1000 = as.integer(round(score * 1000))
    ) |>
    dplyr::select(
      gene_a, gene_b, score, score_1000, is_nqo1_edge,
      nscore, fscore, pscore, ascore, escore, dscore, tscore,
      string_id_a, string_id_b
    )

  log_msg(sprintf("STRING: %d arestas únicas; %d ligadas diretamente a NQO1.",
                  nrow(df), sum(df$is_nqo1_edge)))
  df
}

montar_grafo <- function(arestas) {
  igraph::graph_from_data_frame(
    arestas[, c("gene_a", "gene_b", "score", "score_1000", "is_nqo1_edge")],
    directed = FALSE
  )
}

# ------------------------------------------------------------------------------
# 4) EXPRESSÃO DIFERENCIAL (limma) — OPCIONAL
# ------------------------------------------------------------------------------
carregar_deg <- function() {
  csv <- file.path(TMP, "deg_full.csv")
  if (!file.exists(csv)) {
    rds <- file.path(ROOT, "data", "processed", "resultado_limma_DCM_vs_CTRL.rds")
    if (!file.exists(rds)) {
      log_msg("Resultado do limma não encontrado; seguindo sem log2FC.")
      return(tibble::tibble(
        gene = character(), log2FoldChange = numeric(), padj = numeric(),
        significant = logical(), direction = character()
      ))
    }
    x <- readRDS(rds)
    readr::write_csv(x$deg, csv)
  }

  deg <- readr::read_csv(csv, show_col_types = FALSE, progress = FALSE)
  colunas <- c("gene", "log2FoldChange", "padj", "significant", "direction")
  for (coluna in colunas) {
    if (!coluna %in% names(deg)) deg[[coluna]] <- NA
  }
  deg[, colunas]
}

# ------------------------------------------------------------------------------
# 5) KEGG — VIAS E GENES
# ------------------------------------------------------------------------------
parse_kegg_entry <- function(texto) {
  linhas <- strsplit(texto, "\n", fixed = TRUE)[[1]]
  obter_campo <- function(campo) {
    idx <- which(startsWith(linhas, campo))
    if (!length(idx)) return(NA_character_)
    trimws(sub(paste0("^", campo, "\\s+"), "", linhas[idx[1]]))
  }
  list(name = obter_campo("NAME"),
       class = obter_campo("CLASS"),
       entry = obter_campo("ENTRY"))
}

baixar_kegg <- function() {
  log_msg("Consultando KEGG...")

  link_nqo1 <- cached_text(
    "kegg_nqo1_link_pathway.tsv",
    "https://rest.kegg.jp/link/pathway/hsa:1728"
  )
  linhas <- strsplit(link_nqo1, "\n", fixed = TRUE)[[1]]
  linhas <- linhas[nzchar(linhas)]
  pathway_ids <- unique(sub("^.*\tpath:", "", linhas))

  listagem <- cached_text("kegg_list_hsa.tsv", "https://rest.kegg.jp/list/hsa")
  linhas <- strsplit(listagem, "\n", fixed = TRUE)[[1]]
  id2sym <- list()
  for (linha in linhas) {
    campos <- strsplit(linha, "\t", fixed = TRUE)[[1]]
    if (length(campos) < 2) next
    kid <- campos[1]
    desc <- campos[length(campos)]
    if (!grepl(";", desc, fixed = TRUE)) next
    simbolos <- trimws(strsplit(sub(";.*$", "", desc), ",", fixed = TRUE)[[1]])
    simbolos <- simbolos[nzchar(simbolos)]
    if (!length(simbolos)) next
    id2sym[[kid]] <- simbolos[1]
  }

  genes_por_via <- list()
  for (pid in pathway_ids) {
    texto <- cached_text(
      paste0("kegg_link_hsa_", pid, ".tsv"),
      paste0("https://rest.kegg.jp/link/hsa/", pid)
    )
    linhas <- strsplit(texto, "\n", fixed = TRUE)[[1]]
    linhas <- linhas[nzchar(linhas)]
    kids <- sub("^.*\t", "", linhas)
    genes <- unique(unlist(id2sym[kids], use.names = FALSE))
    genes_por_via[[pid]] <- genes[!is.na(genes)]
  }

  linhas_out <- list()
  for (pid in pathway_ids) {
    info <- parse_kegg_entry(cached_text(
      paste0("kegg_get_", pid, ".txt"),
      paste0("https://rest.kegg.jp/get/", pid)
    ))

    nome <- info$name
    if (!is.na(nome)) {
      nome <- sub("\\s+-\\s+Homo sapiens \\(human\\)$", "", nome)
    }
    classe <- info$class
    if (is.na(classe) || !nzchar(classe)) {
      classe <- if (pid == "hsa01100") {
        "Metabolism; Global and overview maps"
      } else if (pid == "hsa01240") {
        "Metabolism; Metabolism of cofactors and vitamins"
      } else {
        "Não classificada"
      }
    }

    linhas_out[[length(linhas_out) + 1]] <- tibble::tibble(
      pathway_id = pid,
      pathway_name = if (is.na(nome) || !nzchar(nome)) pid else nome,
      class = classe,
      class_group = trimws(strsplit(classe, ";", fixed = TRUE)[[1]][1]),
      n_genes_pathway = length(genes_por_via[[pid]]),
      url = paste0("https://www.kegg.jp/pathway/", pid)
    )
  }

  kegg <- dplyr::bind_rows(linhas_out)
  log_msg(sprintf("KEGG: %d vias; classes: %s.",
                  nrow(kegg), paste(sort(unique(kegg$class_group)), collapse = ", ")))
  list(kegg = kegg, genes_por_via = genes_por_via)
}

# ------------------------------------------------------------------------------
# 6) REACTOME — VIAS, HIERARQUIA E PARTICIPANTES
# ------------------------------------------------------------------------------
extrair_cadeia <- function(ancestrais) {
  if (is.null(ancestrais) || length(ancestrais) == 0) return(character())
  primeiro <- ancestrais[[1]]
  caminho <- if (is.list(primeiro) && !is.null(primeiro$displayName)) {
    ancestrais
  } else {
    primeiro
  }
  valores <- vapply(caminho, function(x) {
    if (!is.null(x$displayName)) as.character(x$displayName) else ""
  }, character(1))
  valores[nzchar(valores)]
}

extrair_genes_reactome <- function(participantes) {
  genes <- character()
  for (p in participantes) {
    refs <- p$refEntities
    if (is.null(refs)) next
    for (ref in refs) {
      dn <- ref$displayName
      if (is.null(dn) || !nzchar(dn)) next
      if (grepl("^(UniProt|ENSEMBL):\\S+\\s+", dn)) {
        genes <- c(genes, sub("^(UniProt|ENSEMBL):\\S+\\s+", "", dn))
      }
    }
  }
  unique(genes)
}

baixar_reactome <- function() {
  log_msg("Consultando Reactome...")

  busca <- cached_json(
    "reactome_search_nqo1.json",
    "https://reactome.org/ContentService/search/query",
    query = list(query = GENE, species = "Homo sapiens", cluster = "true")
  )

  db_id <- NULL
  for (resultado in busca$results) {
    for (entrada in resultado$entries) {
      if (identical(entrada$type, "Protein") &&
          identical(entrada$referenceIdentifier, UNIPROT)) {
        db_id <- as.character(entrada$dbId)
        break
      }
    }
    if (!is.null(db_id)) break
  }
  if (is.null(db_id)) {
    stop("Não consegui mapear NQO1 para uma entidade do Reactome.", call. = FALSE)
  }

  vias_baixas <- cached_json(
    "reactome_nqo1_low_pathways.json",
    paste0("https://reactome.org/ContentService/data/pathways/low/entity/", db_id, "/allForms"),
    query = list(species = SPECIES)
  )

  genes_por_via <- list()
  linhas <- list()
  hierarquia <- list()

  for (via in vias_baixas) {
    vid <- as.character(via$dbId)
    stid <- via$stId
    nome <- via$displayName
    if (is.null(nome) || !nzchar(nome)) nome <- via$name[[1]]

    ancestrais <- cached_json(
      paste0("reactome_ancestors_", vid, ".json"),
      paste0("https://reactome.org/ContentService/data/event/", vid, "/ancestors")
    )
    cadeia <- extrair_cadeia(ancestrais)
    if (!length(cadeia)) cadeia <- nome
    top_level <- cadeia[length(cadeia)]

    participantes <- cached_json(
      paste0("reactome_participants_", vid, ".json"),
      paste0("https://reactome.org/ContentService/data/participants/", vid)
    )
    genes <- extrair_genes_reactome(participantes)
    genes_por_via[[stid]] <- genes

    linhas[[length(linhas) + 1]] <- tibble::tibble(
      pathway_id = stid,
      pathway_name = nome,
      top_level = top_level,
      ancestor_chain = paste(cadeia, collapse = " > "),
      n_participants = length(genes),
      n_ancestors = length(cadeia),
      url = paste0("https://reactome.org/content/detail/", stid)
    )

    for (ordem in seq_along(cadeia)) {
      ordem_evento <- ordem - 1L
      evento_atual <- cadeia[[ordem]]
      direto_atual <- ordem == 1L
      hierarquia[[length(hierarquia) + 1]] <- tibble::tibble(
        pathway_id = stid,
        pathway_name = nome,
        ordem = ordem_evento,
        evento = evento_atual,
        direto = direto_atual,
        top_level = top_level,
        url = paste0("https://reactome.org/content/detail/", stid)
      )
    }
  }

  reactome <- dplyr::bind_rows(linhas)
  hierarquia_df <- dplyr::bind_rows(hierarquia)
  log_msg(sprintf("Reactome: %d vias de baixo nível; top-level: %s.",
                  nrow(reactome),
                  paste(sort(unique(reactome$top_level)), collapse = ", ")))
  list(reactome = reactome, hierarquia = hierarquia_df,
       genes_por_via = genes_por_via, db_id = db_id)
}

# ------------------------------------------------------------------------------
# 7) CRUZAMENTOS E TABELAS
# ------------------------------------------------------------------------------
anotar_nos <- function(g, arestas, kegg, genes_kegg, reactome, genes_reactome, deg) {
  genes_rede <- igraph::V(g)$name

  diretos <- arestas[arestas$is_nqo1_edge, , drop = FALSE]
  score_direto <- setNames(rep(NA_real_, length(genes_rede)), genes_rede)
  for (i in seq_len(nrow(diretos))) {
    a <- diretos$gene_a[i]
    b <- diretos$gene_b[i]
    parceiro <- if (a == GENE) b else a
    score_direto[parceiro] <- diretos$score[i]
  }
  if (nrow(diretos) > 0) score_direto[GENE] <- max(diretos$score, na.rm = TRUE)

  string_ids <- setNames(rep(NA_character_, length(genes_rede)), genes_rede)
  for (i in seq_len(nrow(arestas))) {
    string_ids[arestas$gene_a[i]] <- arestas$string_id_a[i]
    string_ids[arestas$gene_b[i]] <- arestas$string_id_b[i]
  }

  grau <- igraph::degree(g)
  peso <- igraph::strength(g, weights = igraph::E(g)$score)
  intermed <- igraph::betweenness(g, directed = FALSE, weights = NA, normalized = TRUE)
  proximidade <- igraph::closeness(g, normalized = TRUE)
  agrupamento <- igraph::transitivity(g, type = "local")

  linhas <- vector("list", length(genes_rede))
  for (i in seq_along(genes_rede)) {
    no <- genes_rede[i]
    alvo <- alias_gene(no)
    vias_kegg <- names(genes_kegg)[vapply(
      genes_kegg, function(genes) alvo %in% genes, logical(1)
    )]
    vias_reactome <- names(genes_reactome)[vapply(
      genes_reactome, function(genes) alvo %in% genes, logical(1)
    )]

    linhas[[i]] <- tibble::tibble(
      gene = no,
      is_nqo1 = no == GENE,
      string_id = unname(string_ids[no]),
      degree = as.integer(grau[no]),
      weighted_degree = round(unname(peso[no]), 3),
      direct_score_nqo1 = unname(score_direto[no]),
      betweenness = round(unname(intermed[no]), 4),
      closeness = round(unname(proximidade[no]), 4),
      clustering = round(unname(agrupamento[no]), 4),
      n_kegg = length(vias_kegg),
      n_reactome = length(vias_reactome),
      n_pathways_total = length(vias_kegg) + length(vias_reactome),
      kegg_pathway_ids = paste(sort(vias_kegg), collapse = ";"),
      reactome_pathway_ids = paste(sort(vias_reactome), collapse = ";")
    )
  }
  nodes <- dplyr::bind_rows(linhas)

  if (nrow(deg) > 0) {
    nodes <- dplyr::left_join(nodes, deg, by = "gene")
  } else {
    nodes$log2FoldChange <- NA_real_
    nodes$padj <- NA_real_
    nodes$significant <- NA
    nodes$direction <- NA_character_
  }

  kegg_nome <- setNames(kegg$pathway_name, kegg$pathway_id)
  react_nome <- setNames(reactome$pathway_name, reactome$pathway_id)

  nomes_vias <- function(ids, dicionario) {
    if (is.na(ids) || !nzchar(ids)) return("")
    valores <- strsplit(ids, ";", fixed = TRUE)[[1]]
    paste(vapply(valores, function(x) {
      if (x %in% names(dicionario)) unname(dicionario[x]) else x
    }, character(1)), collapse = "; ")
  }

  nodes$kegg_pathway_names <- vapply(
    nodes$kegg_pathway_ids, nomes_vias, character(1), dicionario = kegg_nome
  )
  nodes$reactome_pathway_names <- vapply(
    nodes$reactome_pathway_ids, nomes_vias, character(1), dicionario = react_nome
  )
  nodes$medido <- !is.na(nodes$log2FoldChange)

  nodes |>
    dplyr::arrange(dplyr::desc(is_nqo1), dplyr::desc(degree)) |>
    dplyr::mutate(gene = as.character(gene)) |>
    as.data.frame(stringsAsFactors = FALSE)
}

completar_kegg_reactome <- function(tab, genes_por_via, parceiros) {
  tab$string_partners <- vapply(tab$pathway_id, function(pid) {
    genes <- genes_por_via[[pid]]
    if (is.null(genes)) return("")
    presentes <- parceiros[vapply(
      parceiros, function(p) alias_gene(p) %in% genes, logical(1)
    )]
    paste(sort(presentes), collapse = ";")
  }, character(1))

  tab$n_string_partners <- ifelse(
    tab$string_partners == "", 0L,
    lengths(strsplit(tab$string_partners, ";", fixed = TRUE))
  )

  tab |>
    dplyr::arrange(dplyr::desc(n_string_partners), pathway_name) |>
    as.data.frame(stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------------------
# 8) FIGURAS
# ------------------------------------------------------------------------------
posicoes_radiais <- function(g, nodes) {
  parceiros <- igraph::V(g)$name[igraph::V(g)$name != GENE]
  ordem <- nodes$gene[nodes$gene != GENE]
  ordem <- ordem[order(-nodes$degree[nodes$gene != GENE], nodes$gene[nodes$gene != GENE])]
  parceiros <- c(intersect(ordem, parceiros), setdiff(parceiros, ordem))
  angulos <- seq(pi / 2, pi / 2 - 2 * pi, length.out = length(parceiros) + 1)
  angulos <- angulos[seq_len(length(parceiros))]
  tibble::tibble(
    gene = c(GENE, parceiros),
    x = c(0, cos(angulos)),
    y = c(0, sin(angulos))
  )
}

figura_rede <- function(g, nodes, arestas) {
  log_msg("Desenhando a rede PPI...")
  coords <- posicoes_radiais(g, nodes) |>
    dplyr::left_join(nodes, by = "gene")

  coords_from <- coords[, c("gene", "x", "y")]
  coords_to <- coords[, c("gene", "x", "y")]
  names(coords_to) <- c("gene", "xend", "yend")

  edges <- igraph::as_data_frame(g, what = "edges") |>
    dplyr::left_join(coords_from, by = c("from" = "gene")) |>
    dplyr::left_join(coords_to, by = c("to" = "gene")) |>
    dplyr::mutate(
      is_nqo1_edge = from == GENE | to == GENE,
      edge_width = ifelse(
        is_nqo1_edge,
        0.6 + 2.2 * pmax(0, (score - 0.70) / 0.30),
        0.45
      )
    )

  lim_fc <- suppressWarnings(max(abs(coords$log2FoldChange), na.rm = TRUE))
  if (!is.finite(lim_fc)) lim_fc <- 1
  lim_fc <- max(1, ceiling(lim_fc * 10) / 10)

  top_hubs <- nodes |>
    dplyr::filter(gene != GENE) |>
    dplyr::arrange(dplyr::desc(degree)) |>
    dplyr::slice_head(n = 8) |>
    dplyr::pull(gene)

  coords <- coords |>
    dplyr::mutate(
      top_hub = gene %in% top_hubs,
      label_hjust = ifelse(x > 0.05, 0, ifelse(x < -0.05, 1, 0.5)),
      label_x = x * 1.105,
      label_y = y * 1.105
    )

  p <- ggplot() +
    geom_segment(
      data = dplyr::filter(edges, !is_nqo1_edge),
      aes(x = x, y = y, xend = xend, yend = yend),
      color = "grey55", alpha = 0.10, linewidth = 0.25
    ) +
    geom_segment(
      data = dplyr::filter(edges, is_nqo1_edge),
      aes(x = x, y = y, xend = xend, yend = yend, linewidth = edge_width),
      color = "grey25", alpha = 0.55
    ) +
    scale_linewidth_identity(guide = "none") +
    geom_point(
      data = dplyr::filter(coords, gene != GENE, is.na(log2FoldChange)),
      aes(x = x, y = y, size = degree),
      shape = 21, fill = "#c9c9c9", color = "#777777", stroke = 0.4
    ) +
    geom_point(
      data = dplyr::filter(coords, gene != GENE, !is.na(log2FoldChange)),
      aes(x = x, y = y, size = degree, fill = log2FoldChange),
      shape = 21, color = "#222222", stroke = 0.5
    ) +
    geom_point(
      data = dplyr::filter(coords, gene == GENE),
      aes(x = x, y = y, fill = log2FoldChange),
      shape = 21, size = 16, color = "black", stroke = 1.5
    ) +
    geom_text(
      data = dplyr::filter(coords, gene != GENE),
      aes(x = label_x, y = label_y, label = gene,
          hjust = label_hjust, fontface = ifelse(top_hub, "bold", "plain")),
      size = 3, color = "#111111"
    ) +
    geom_text(
      data = dplyr::filter(coords, gene == GENE),
      aes(x = x, y = y, label = gene),
      color = "white", fontface = "bold", size = 4
    ) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-lim_fc, lim_fc),
      na.value = "#c9c9c9", name = "log2FC\n(DCM vs CTRL)"
    ) +
    scale_size_continuous(
      range = c(2, 10), breaks = c(5, 15, 25),
      name = "Grau"
    ) +
    coord_equal(clip = "off", xlim = c(-1.32, 1.32), ylim = c(-1.25, 1.25)) +
    labs(
      title = "Rede PPI centrada em NQO1 no STRING",
      subtitle = sprintf(
        "Homo sapiens | score ≥ %d | %d nós e %d arestas | %d parceiros diretos",
        STRING_SCORE, igraph::gorder(g), igraph::gsize(g), igraph::degree(g)[GENE]
      ),
      caption = paste0(
        "Anel externo: parceiros diretos de NQO1. Cor: log2FC na DCM ",
        "(azul = down, vermelho = up; cinza = gene não medido no dataset).\n",
        "Arestas escuras: interação NQO1–parceiro; arestas claras: interação entre parceiros. ",
        "Tamanho do nó: grau na rede."
      )
    ) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 10.5, hjust = 0.5, color = "#444444"),
      plot.caption = element_text(size = 8.2, color = "#555555", hjust = 0)
    )

  ggsave(file.path(FIG, "nqo1_rede_ppi_string.png"), p,
         width = 16, height = 14, dpi = 300, bg = "white")
  log_msg("Figura da rede salva.")
}

quebrar_texto <- function(texto, largura = 42) {
  stringr::str_wrap(texto, width = largura)
}

figura_kegg <- function(kegg) {
  log_msg("Desenhando barras KEGG...")
  df <- kegg |>
    dplyr::arrange(n_string_partners, pathway_name) |>
    dplyr::mutate(
      rotulo = paste0(quebrar_texto(pathway_name, 48), "\n[", pathway_id, "]"),
      classe_plot = factor(class_group,
                           levels = c("Metabolism", "Human Diseases", "Não classificada"))
    )

  cores <- c("Metabolism" = "#2A9D8F",
             "Human Diseases" = "#E76F51",
             "Não classificada" = "#888888")

  p <- ggplot(df, aes(x = n_string_partners, y = reorder(rotulo, n_string_partners),
                      fill = classe_plot)) +
    geom_col(width = 0.68, color = "#333333", linewidth = 0.3) +
    geom_text(
      aes(label = sprintf("%d parceiros  (%d genes na via)",
                          n_string_partners, n_genes_pathway)),
      hjust = 0, nudge_x = 0.15, size = 3, color = "#333333"
    ) +
    scale_fill_manual(values = cores, drop = FALSE, name = "Classe KEGG") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "Vias KEGG que contêm NQO1",
      subtitle = paste0(
        "Barras = quantos parceiros da rede STRING também estão anotados na via. ",
        "Entre parênteses: total de genes da via no KEGG."
      ),
      x = "Nº de parceiros diretos de NQO1 (STRING) presentes na via",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(size = 9, color = "#555555"),
      axis.text.y = element_text(size = 8.3),
      axis.text.x = element_text(size = 9)
    )

  ggsave(file.path(FIG, "nqo1_kegg_pathways.png"), p,
         width = 13, height = 7.2, dpi = 300, bg = "white")
  log_msg("Figura KEGG salva.")
}

figura_reactome <- function(reactome) {
  log_msg("Desenhando barras Reactome...")
  df <- reactome |>
    dplyr::arrange(n_string_partners, pathway_name) |>
    dplyr::mutate(
      rotulo = paste0(quebrar_texto(pathway_name, 44), "\n[", pathway_id, "]")
    )

  palette <- c(
    "Metabolism" = "#264653",
    "Chromatin organization" = "#8E7DBE",
    "Cellular responses to stimuli" = "#E9C46A",
    "Cellular responses to stress" = "#F4A261",
    "Cellular response to chemical stress" = "#E76F51"
  )

  p <- ggplot(df, aes(x = n_string_partners, y = reorder(rotulo, n_string_partners),
                      fill = top_level)) +
    geom_col(width = 0.6, color = "#333333", linewidth = 0.3) +
    geom_text(
      aes(label = sprintf("%d parceiros  (%d genes/participantes)",
                          n_string_partners, n_participants)),
      hjust = 0, nudge_x = 0.12, size = 3, color = "#333333"
    ) +
    scale_fill_manual(values = palette, na.value = "#888888",
                      name = "Top-level Reactome") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
      title = "Vias Reactome que contêm NQO1",
      subtitle = paste0(
        "Barras = parceiros STRING também presentes na via. ",
        "Cor = via de nível superior (top-level) do Reactome."
      ),
      x = "Nº de parceiros diretos de NQO1 (STRING) também anotados na via",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(size = 9, color = "#555555"),
      axis.text.y = element_text(size = 8.5),
      axis.text.x = element_text(size = 9)
    )

  ggsave(file.path(FIG, "nqo1_reactome_pathways.png"), p,
         width = 13, height = 5.4, dpi = 300, bg = "white")
  log_msg("Figura Reactome salva.")
}

figura_heatmap <- function(nodes, kegg, reactome) {
  log_msg("Desenhando heatmap de vias...")
  df <- nodes |>
    dplyr::filter(is_nqo1 | n_pathways_total > 0) |>
    dplyr::arrange(dplyr::desc(is_nqo1), dplyr::desc(n_pathways_total),
                   dplyr::desc(degree))
  genes <- df$gene

  ids_kegg <- sort(unique(kegg$pathway_id))
  ids_react <- sort(unique(reactome$pathway_id))
  ids <- c(ids_kegg, ids_react)

  rotulos <- c(
    setNames(
      paste0(ids_kegg, "\n", quebrar_texto(
        kegg$pathway_name[match(ids_kegg, kegg$pathway_id)], 26
      )),
      ids_kegg
    ),
    setNames(
      paste0(ids_react, "\n", quebrar_texto(
        reactome$pathway_name[match(ids_react, reactome$pathway_id)], 26
      )),
      ids_react
    )
  )

  heat <- list()
  for (gene in genes) {
    linha <- df[df$gene == gene, , drop = FALSE]
    k_ids <- strsplit(linha$kegg_pathway_ids, ";", fixed = TRUE)[[1]]
    r_ids <- strsplit(linha$reactome_pathway_ids, ";", fixed = TRUE)[[1]]
    k_ids <- k_ids[nzchar(k_ids)]
    r_ids <- r_ids[nzchar(r_ids)]

    for (pid in ids_kegg) {
      heat[[length(heat) + 1]] <- tibble::tibble(
        gene = gene, pathway_id = pid, grupo = "KEGG",
        valor = as.integer(pid %in% k_ids)
      )
    }
    for (pid in ids_react) {
      heat[[length(heat) + 1]] <- tibble::tibble(
        gene = gene, pathway_id = pid, grupo = "Reactome",
        valor = as.integer(pid %in% r_ids)
      )
    }
  }
  heat <- dplyr::bind_rows(heat) |>
    dplyr::mutate(
      gene = factor(gene, levels = rev(genes)),
      pathway_id = factor(pathway_id, levels = ids),
      grupo = factor(grupo, levels = c("KEGG", "Reactome"))
    )

  altura <- max(5.5, 0.28 * length(genes) + 2.2)

  p <- ggplot(heat, aes(x = pathway_id, y = gene, fill = factor(valor))) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(data = dplyr::filter(heat, valor == 1),
              label = "•", color = "white", size = 3) +
    facet_grid(~ grupo, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(`0` = "#f7f7f7", `1` = "#1f6f8b"),
                      guide = "none") +
    scale_x_discrete(labels = rotulos) +
    scale_y_discrete(labels = function(x) ifelse(x == GENE, paste0(">> ", x, " <<"), x)) +
    labs(
      title = "Parceiros de NQO1 que compartilham vias KEGG/Reactome",
      x = "Via (ID e nome abreviado)",
      y = "Proteína"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7.2),
      axis.text.y = element_text(size = 8.2),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid = element_blank(),
      panel.spacing.x = unit(0.4, "lines")
    )

  ggsave(file.path(FIG, "nqo1_heatmap_vias_parceiros.png"), p,
         width = 13.5, height = altura, dpi = 300, bg = "white")
  log_msg("Heatmap salvo.")
}

figura_kegg_diagrama <- function() {
  log_msg("Baixando e destacando o mapa KEGG hsa05208...")
  tryCatch({
    bytes <- cached_raw(
      "kegg_hsa05208.png",
      "https://rest.kegg.jp/get/hsa05208/image"
    )
    img <- magick::image_read(bytes)

    html <- cached_text(
      "kegg_hsa05208_highlight.html",
      "https://www.kegg.jp/kegg-bin/show_pathway?hsa05208+1728"
    )
    tags <- regmatches(html, gregexpr("<area\\b[^>]*>", html, perl = TRUE))[[1]]
    tags <- tags[grepl("NQO1", tags, fixed = TRUE) & grepl("coords=", tags, fixed = TRUE)]

    coords <- lapply(tags, function(tag) {
      m <- regmatches(tag, regexpr('coords="([^"]+)"', tag, perl = TRUE))
      if (!length(m)) return(NULL)
      nums <- as.integer(regmatches(m, gregexpr("\\d+", m, perl = TRUE))[[1]])
      if (length(nums) >= 4) nums[1:4] else NULL
    })
    coords <- Filter(Negate(is.null), coords)

    destino <- file.path(FIG, "nqo1_kegg_diagram_hsa05208.png")
    if (!length(coords)) {
      log_msg("Não encontrei as coordenadas do NQO1 no HTML do KEGG; salvando mapa sem destaque.")
      magick::image_write(img, destino)
      return(invisible(NULL))
    }

    img <- magick::image_draw(img)
    for (i in seq_along(coords)) {
      r <- coords[[i]]
      graphics::rect(r[1] - 7, r[2] - 7, r[3] + 7, r[4] + 7,
                     border = "#DC1414", lwd = 5)
      if (i == 1) {
        graphics::text(r[1] - 7, max(0, r[2] - 42), "NQO1",
                       col = "#DC1414", cex = 1.6, adj = c(0, 1), font = 2)
      }
    }
    grDevices::dev.off()
    magick::image_write(img, destino)
    log_msg(sprintf("Mapa KEGG salvo com %d destaque(s) para NQO1.", length(coords)))
  }, error = function(e) {
    log_msg("Não foi possível gerar o mapa KEGG: ", conditionMessage(e))
  })
}

figura_reactome_diagrama <- function(reactome) {
  log_msg("Baixando o diagrama Reactome da via NFE2L2...")
  tryCatch({
    if (!"R-HSA-9818027" %in% reactome$pathway_id) return(invisible(NULL))
    bytes <- cached_raw(
      "reactome_R-HSA-9818027.png",
      "https://reactome.org/ContentService/exporter/diagram/R-HSA-9818027.png"
    )
    writeBin(bytes, file.path(FIG, "nqo1_reactome_diagram_nfe2l2.png"))
    log_msg("Diagrama Reactome salvo.")
  }, error = function(e) {
    log_msg("Não foi possível baixar o diagrama Reactome: ", conditionMessage(e))
  })
}

# ------------------------------------------------------------------------------
# 9) MÉTRICAS
# ------------------------------------------------------------------------------
salvar_metricas <- function(nodes, kegg, reactome, g) {
  log_msg("Salvando métricas do resumo...")
  nqo <- nodes[nodes$is_nqo1, , drop = FALSE][1, , drop = FALSE]

  metricas <- tibble::tibble(
    Métrica = c(
      "Gene alvo",
      "UniProt",
      "Espécie",
      "Rede STRING — nós",
      "Rede STRING — arestas",
      "Parceiros diretos de NQO1",
      "Score mínimo STRING",
      "Vias KEGG com NQO1",
      "Vias Reactome (baixo nível) com NQO1",
      "NQO1 log2FC (DCM vs CTRL)",
      "NQO1 padj",
      "NQO1 direção"
    ),
    Valor = c(
      GENE,
      UNIPROT,
      "Homo sapiens (9606)",
      as.character(igraph::gorder(g)),
      as.character(igraph::gsize(g)),
      as.character(igraph::degree(g)[GENE]),
      as.character(STRING_SCORE),
      as.character(nrow(kegg)),
      as.character(nrow(reactome)),
      sprintf("%.3f", nqo$log2FoldChange),
      sprintf("%.2e", nqo$padj),
      as.character(nqo$direction)
    )
  )

  readr::write_csv(metricas, file.path(TAB, "nqo1_resumo_metricas.csv"))
  log_msg("Métricas salvas em output/tables/nqo1_resumo_metricas.csv.")
  log_msg("A interpretação está integrada em RELATORIO_CIENTIFICO.md, seção 4.")
}

# ------------------------------------------------------------------------------
# 10) MAIN
# ------------------------------------------------------------------------------
main <- function() {
  log_msg("Iniciando análise NQO1: STRING + KEGG + Reactome (R).")

  arestas <- baixar_rede_string()
  g <- montar_grafo(arestas)
  deg <- carregar_deg()

  kegg_info <- baixar_kegg()
  kegg <- kegg_info$kegg
  genes_kegg <- kegg_info$genes_por_via

  reactome_info <- baixar_reactome()
  reactome <- reactome_info$reactome
  hierarquia <- reactome_info$hierarquia
  genes_reactome <- reactome_info$genes_por_via

  parceiros <- setdiff(igraph::V(g)$name, GENE)
  kegg <- completar_kegg_reactome(kegg, genes_kegg, parceiros)
  reactome <- completar_kegg_reactome(reactome, genes_reactome, parceiros)

  nodes <- anotar_nos(
    g, arestas, kegg, genes_kegg, reactome, genes_reactome, deg
  )

  # Tabelas.
  readr::write_csv(arestas, file.path(TAB, "nqo1_string_ppi_arestas.csv"), na = "")
  readr::write_csv(nodes, file.path(TAB, "nqo1_string_ppi_nos.csv"), na = "")
  readr::write_csv(kegg, file.path(TAB, "nqo1_kegg_pathways.csv"), na = "")
  readr::write_csv(reactome, file.path(TAB, "nqo1_reactome_pathways.csv"), na = "")
  readr::write_csv(hierarquia, file.path(TAB, "nqo1_reactome_hierarquia.csv"), na = "")
  log_msg("Tabelas CSV salvas em output/tables/.")

  # Figuras.
  figura_rede(g, nodes, arestas)
  figura_kegg(kegg)
  figura_reactome(reactome)
  figura_heatmap(nodes, kegg, reactome)
  figura_kegg_diagrama()
  figura_reactome_diagrama(reactome)

  # Métricas.
  salvar_metricas(nodes, kegg, reactome, g)

  log_msg("Concluído.")
  log_msg(sprintf("Rede: %d nós, %d arestas, %d parceiros diretos de NQO1.",
                  igraph::gorder(g), igraph::gsize(g), igraph::degree(g)[GENE]))
  log_msg(sprintf("KEGG: %d vias | Reactome: %d vias.", nrow(kegg), nrow(reactome)))
}

main()
