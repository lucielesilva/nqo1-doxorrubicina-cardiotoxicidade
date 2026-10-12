#!/usr/bin/env Rscript

# ==============================================================================
# Análise integrativa in silico da regulação do NQO1 e da cardiotoxicidade
# da doxorrubicina na cardiomiopatia dilatada
# Pipeline analítico único
# ==============================================================================
# Este arquivo reúne todo o código do projeto, em três partes:
#
#   PARTE A. Análise upstream (PASSO 1 + PASSO 2)
#     Expressão diferencial do GSE116250 com limma, caracterização da via AMPK
#     (KEGG hsa04152) por GSEA e GSVA e seleção do gene NQO1.
#
#   PARTE B. Análise downstream (PASSO 3)
#     Leitura do DGB (CMAP, LINCS L1000, CREEDS), ranqueamento dos fármacos que
#     modulam o NQO1, figuras, PubChem e rede PPI (STRING).
#
#   PARTE C. Auditoria de dados
#     Integridade, valores ausentes, duplicatas, outliers e consistência entre
#     amostras e metadados.
#
# Estrutura de pastas:
#   data/raw/                         dados brutos (GEO e STRING)
#   output/rdata/                     matriz pré-processada
#   output/tables/{ampk,nqo1,pubchem,ppi}/
#   output/figures/{ampk,nqo1_drug_regulation,pubchem,ppi}/
#   output/audit/                     logs e relatórios de auditoria
#   output/archive/                   resultados antigos preservados
#
# Como rodar no PowerShell, a partir da raiz do projeto:
#   & 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' analise_nqo1_cardiotoxicidade.R
#
# Dependências obrigatórias: dplyr, ggplot2, stringr, openxlsx, igraph, tidyr,
# data.table, readr, limma.
# Dependências opcionais: fgsea (GSEA) e GSVA (escore por amostra).
#
# Consulta PubChem:
#   Por padrão, o script usa os resultados PubChem já salvos em
#   output/tables/pubchem. Para refazer a consulta externa, defina
#   RUN_PUBCHEM_QUERY e NCBI_API_KEY no ambiente. A chave nunca fica no código.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(42)

# ------------------------------------------------------------------------------
# 0. Caminhos e utilitarios
# ------------------------------------------------------------------------------

get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT <- get_project_root()
setwd(ROOT)

DIRS <- list(
  data_raw = file.path(ROOT, "data", "raw"),
  rdata = file.path(ROOT, "output", "rdata"),
  tables_ampk = file.path(ROOT, "output", "tables", "ampk"),
  tables_pubchem = file.path(ROOT, "output", "tables", "pubchem"),
  tables_nqo1 = file.path(ROOT, "output", "tables", "nqo1"),
  tables_ppi = file.path(ROOT, "output", "tables", "ppi"),
  figures_ampk = file.path(ROOT, "output", "figures", "ampk"),
  figures_pubchem = file.path(ROOT, "output", "figures", "pubchem"),
  figures_nqo1 = file.path(ROOT, "output", "figures", "nqo1_drug_regulation"),
  figures_ppi = file.path(ROOT, "output", "figures", "ppi"),
  audit = file.path(ROOT, "output", "audit"),
  logs = file.path(ROOT, "output", "audit", "logs")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

AUDIT_EVENTS <- data.frame(
  time = character(),
  level = character(),
  step = character(),
  message = character()
)

audit_log <- function(level, step, message) {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %-5s | %-18s | %s", stamp, level, step, message)
  message(line)
  log_file <- file.path(DIRS$logs, sprintf("execucao_%s.log", format(Sys.Date(), "%Y%m%d")))
  cat(line, "\n", file = log_file, append = TRUE)
  AUDIT_EVENTS <<- rbind(
    AUDIT_EVENTS,
    data.frame(time = stamp, level = level, step = step, message = message)
  )
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Pacote obrigatório ausente: ", pkg, ". Instale antes de rodar este script.")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

require_pkg("dplyr")
require_pkg("ggplot2")
require_pkg("stringr")
require_pkg("openxlsx")
require_pkg("igraph")
require_pkg("tidyr")
require_pkg("data.table")
require_pkg("readr")
require_pkg("limma")

`%+%` <- function(a, b) paste0(a, b)

find_first_existing <- function(paths, required = TRUE, label = "arquivo") {
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0) {
    normalizePath(existing[1], winslash = "/", mustWork = TRUE)
  } else if (required) {
    stop("Não encontrei ", label, ". Caminhos testados: ", paste(paths, collapse = " | "))
  } else {
    NA_character_
  }
}

save_plot <- function(plot, filename, width = 9, height = 6) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  audit_log("INFO", "plot", paste("Figura salva:", normalizePath(filename, winslash = "/", mustWork = FALSE)))
}

read_csv_flexible <- function(path) {
  audit_log("INFO", "read", paste("Lendo", path))
  df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(df) <- trimws(names(df))
  df
}

standardize_names <- function(x) {
  x |>
    stringr::str_replace_all("\\s+", "_") |>
    stringr::str_replace_all("-", "_") |>
    stringr::str_replace_all("\\.", "_")
}

decimal_to_numeric <- function(x) {
  x <- as.character(x)
  has_comma <- stringr::str_detect(x, ",")
  x <- ifelse(
    has_comma,
    stringr::str_replace(stringr::str_replace_all(x, "\\.", ""), ",", "."),
    x
  )
  suppressWarnings(as.numeric(x))
}

q_to_neglog10 <- function(q_value, cap = 50) {
  y <- -log10(pmax(q_value, 10^-cap, na.rm = TRUE))
  pmin(y, cap)
}

extract_dose_um <- function(dose) {
  dose <- as.character(dose)
  dose <- stringr::str_replace_all(dose, ",", ".")
  value <- suppressWarnings(as.numeric(stringr::str_extract(dose, "\\d+(?:\\.\\d+)?")))
  unit <- stringr::str_to_lower(stringr::str_extract(dose, "nm|µm|μm|um|mm"))
  dplyr::case_when(
    is.na(value) ~ NA_real_,
    unit == "nm" ~ value / 1000,
    unit == "mm" ~ value * 1000,
    TRUE ~ value
  )
}

# ==============================================================================
# PARTE A. Análise upstream (PASSO 1 + PASSO 2)
# ==============================================================================
# PASSO 1: análise exploratória da via AMPK (KEGG hsa04152) na DCM (GSE116250).
# PASSO 2: expressão diferencial do transcriptoma completo, que identifica o
# NQO1 como gene downregulado (log2FC aproximadamente -1,69).
#
# Os dados pre-processados ficam em output/rdata/preprocessed_microarray.RData,
# com a matriz log2 `expr_matrix_norm` (41842 genes x 51 amostras) e o objeto
# `pheno_data` (DCM = 37, CTRL = 14). O contraste DCM vs CTRL exclui a ICM e
# ajusta por sexo e idade.
# ==============================================================================

run_upstream_analysis <- function() {
  audit_log("INFO", "upstream", "Início da análise upstream (PASSO 1 + PASSO 2).")

  rdata_path <- file.path(DIRS$rdata, "preprocessed_microarray.RData")
  if (!file.exists(rdata_path)) {
    stop("Arquivo pre-processado ausente: ", rdata_path)
  }
  env <- new.env()
  load(rdata_path, envir = env)
  expr_norm <- env$expr_matrix_norm
  pheno <- env$pheno_data

  if (anyDuplicated(rownames(expr_norm))) {
    rownames(expr_norm) <- make.unique(rownames(expr_norm))
  }
  audit_log("INFO", "upstream", "Matriz pre-processada: " %+% nrow(expr_norm) %+%
              " genes x " %+% ncol(expr_norm) %+% " amostras")

  pheno <- pheno |>
    dplyr::mutate(
      group = factor(dplyr::case_when(
        `disease:ch1` == "non-failing" ~ "CTRL",
        `disease:ch1` == "dilated cardiomyopathy" ~ "DCM",
        TRUE ~ NA_character_
      ), levels = c("CTRL", "DCM")),
      sex = `Sex:ch1`,
      age = as.numeric(`age:ch1`)
    ) |>
    dplyr::filter(title %in% colnames(expr_norm))
  pheno <- pheno[match(colnames(expr_norm), pheno$title), ]

  audit_log("INFO", "upstream", "Amostras: CTRL = " %+% sum(pheno$group == "CTRL") %+%
              " | DCM = " %+% sum(pheno$group == "DCM") %+%
              " (ICM excluida do contraste)")

  # PASSO 2: expressão diferencial do transcriptoma completo.
  audit_log("INFO", "upstream", "Design: ~ 0 + group + sex + age")
  design <- model.matrix(~ 0 + group + sex + age, data = pheno)
  colnames(design) <- sub("^group", "", colnames(design))

  fit <- limma::lmFit(expr_norm, design)
  fit <- limma::contrasts.fit(fit, limma::makeContrasts("DCM-CTRL", levels = design))
  fit <- limma::eBayes(fit, trend = TRUE)

  tt <- limma::topTable(fit, number = Inf, sort.by = "none")
  deg_full <- data.frame(
    gene = rownames(tt),
    log2FoldChange = tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val,
    B = tt$B,
    stringsAsFactors = FALSE,
    row.names = NULL
  ) |>
    dplyr::mutate(
      neg_log10_padj = -log10(pmax(padj, 1e-300)),
      direction = ifelse(log2FoldChange > 0, "up", "down"),
      significant = padj < 0.05
    ) |>
    dplyr::select(gene, log2FoldChange, AveExpr, t, pvalue, padj, B,
                  direction, significant, neg_log10_padj)

  readr::write_csv(deg_full, file.path(DIRS$tables_ampk, "expressao_diferencial_completa.csv"))
  audit_log("INFO", "upstream", "expressao_diferencial_completa.csv: " %+% nrow(deg_full) %+% " genes")

  nqo1 <- deg_full |> dplyr::filter(gene == "NQO1")
  if (nrow(nqo1) == 0) {
    audit_log("WARN", "upstream", "NQO1 não encontrado no DEG.")
  } else {
    audit_log("INFO", "upstream", sprintf("PASSO 2 | NQO1: log2FC = %.4f | padj = %.3e | %s",
                                          nqo1$log2FoldChange[1], nqo1$padj[1], nqo1$direction[1]))
  }
  readr::write_csv(nqo1, file.path(DIRS$tables_ampk, "expressao_nqo1_diferencial.csv"))

  # PASSO 1: via AMPK (hsa04152, KEGG).
  ampk_path <- file.path(DIRS$tables_ampk, "genes_da_via_ampk_kegg.csv")
  if (!file.exists(ampk_path)) {
    audit_log("WARN", "upstream", "genes_da_via_ampk_kegg.csv ausente; via AMPK ignorada.")
    return(invisible(NULL))
  }
  ampk_genes <- readr::read_csv(ampk_path, show_col_types = FALSE)
  ampk_symbols <- unique(ampk_genes$symbol)
  audit_log("INFO", "upstream", "PASSO 1 | Via AMPK (hsa04152): " %+%
              length(ampk_symbols) %+% " genes no KEGG")

  deg_ampk <- deg_full |>
    dplyr::filter(gene %in% ampk_symbols) |>
    dplyr::left_join(
      ampk_genes |> dplyr::select(symbol, description) |> dplyr::distinct(symbol, .keep_all = TRUE),
      by = c("gene" = "symbol")
    )
  readr::write_csv(deg_ampk, file.path(DIRS$tables_ampk, "expressao_diferencial_via_ampk.csv"))
  audit_log("INFO", "upstream", "expressao_diferencial_via_ampk.csv: " %+% nrow(deg_ampk) %+%
              " genes da via AMPK mapeados")

  # GSEA da via AMPK (fgsea).
  if (requireNamespace("fgsea", quietly = TRUE)) {
    ranked <- deg_full |>
      dplyr::filter(!is.na(log2FoldChange)) |>
      dplyr::arrange(dplyr::desc(log2FoldChange))
    stats <- stats::setNames(ranked$log2FoldChange, ranked$gene)

    gsea_ampk <- fgsea::fgseaSimple(
      pathways = list(AMPK_hsa04152 = ampk_symbols),
      stats = stats, nperm = 10000
    )
    gsea_out <- data.frame(
      pathway = "KEGG_AMPK_SIGNALING_PATHWAY_HSA04152",
      pval = gsea_ampk$pval,
      padj = gsea_ampk$padj,
      log2err = NA_real_,
      ES = gsea_ampk$ES,
      NES = gsea_ampk$NES,
      size = gsea_ampk$size,
      leadingEdge = vapply(gsea_ampk$leadingEdge, paste, collapse = "; ", FUN.VALUE = character(1))
    )
    readr::write_csv(gsea_out, file.path(DIRS$tables_ampk, "gsea_via_ampk.csv"))
    audit_log("INFO", "upstream", sprintf("PASSO 1 | GSEA AMPK: NES = %.3f | padj = %.3f",
                                          gsea_out$NES[1], gsea_out$padj[1]))
  } else {
    audit_log("WARN", "upstream", "fgsea indisponivel; gsea_via_ampk.csv não gerado.")
  }

  # GSVA da via AMPK por amostra.
  if (requireNamespace("GSVA", quietly = TRUE)) {
    param <- GSVA::gsvaParam(
      exprData = expr_norm,
      geneSets = list(AMPK = ampk_symbols),
      kcdf = "Gaussian"
    )
    gsva_obj <- GSVA::gsva(param)
    scores <- gsva_obj
    gsva_df <- data.frame(
      sample = colnames(scores),
      AMPK_score = as.numeric(scores["AMPK", ]),
      group = pheno$group[match(colnames(scores), pheno$title)]
    )
    readr::write_csv(gsva_df, file.path(DIRS$tables_ampk, "gsva_via_ampk.csv"))
    audit_log("INFO", "upstream", sprintf("PASSO 1 | GSVA AMPK: %d amostras (DCM = %d, CTRL = %d)",
                                          nrow(gsva_df), sum(gsva_df$group == "DCM"),
                                          sum(gsva_df$group == "CTRL")))
  } else {
    audit_log("WARN", "upstream", "GSVA indisponivel; gsva_via_ampk.csv não gerado.")
  }

  audit_log("INFO", "upstream", "Análise upstream concluída.")
  invisible(NULL)
}

# ==============================================================================
# PARTE B. Análise downstream (PASSO 3)
# ==============================================================================
# Leitura do DGB, ranqueamento dos fármacos moduladores de NQO1, figuras,
# PubChem e rede PPI.
# ==============================================================================

read_nqo1_xlsx <- function(path) {
  sheets <- openxlsx::getSheetNames(path)
  audit_log("INFO", "xlsx", paste("Abas encontradas:", paste(sheets, collapse = ", ")))

  out <- lapply(sheets, function(sheet) {
    df <- openxlsx::read.xlsx(path, sheet = sheet)
    if (ncol(df) > 0 && grepl("^X\\d+$", names(df)[1])) {
      df <- df[, -1, drop = FALSE]
    }
    names(df) <- standardize_names(names(df))

    if (!"Drug_Name" %in% names(df)) stop("Aba sem coluna Drug Name: ", sheet)
    if (!"p_value" %in% names(df)) stop("Aba sem coluna p-value: ", sheet)
    if (!"q_value" %in% names(df)) stop("Aba sem coluna q-value: ", sheet)
    if (!"Fold_Change" %in% names(df)) stop("Aba sem coluna Fold Change: ", sheet)

    source <- dplyr::case_when(
      stringr::str_detect(sheet, "cmap") ~ "CMAP",
      stringr::str_detect(sheet, "l1000") ~ "LINCS L1000",
      stringr::str_detect(sheet, "creeds") ~ "CREEDS",
      TRUE ~ sheet
    )
    expected_direction <- ifelse(stringr::str_detect(sheet, "_up$"), "up", "down")

    df |>
      dplyr::mutate(
        dplyr::across(dplyr::any_of(c("p_value", "q_value", "Fold_Change")), as.character),
        source_database = source,
        source_sheet = sheet,
        expected_direction = expected_direction
      )
  })

  combined <- dplyr::bind_rows(out)

  combined |>
    dplyr::filter(!is.na(Drug_Name), Drug_Name != "") |>
    dplyr::mutate(
      Drug_Name = stringr::str_squish(as.character(Drug_Name)),
      Drug_Name_Norm = stringr::str_to_lower(stringr::str_squish(as.character(Drug_Name))),
      p_value = decimal_to_numeric(p_value),
      q_value = decimal_to_numeric(q_value),
      Fold_Change = decimal_to_numeric(Fold_Change),
      # Fold_Change do DGB e um escore sinalizado (sinal = direção), não uma
      # razão de fold-change clássica. Mantemos o nome por compatibilidade e
      # adicionamos a semântica explicita (effect_score) e uma versão log2.
      effect_score = Fold_Change,
      log2FC_derived = sign(Fold_Change) * log2(abs(Fold_Change) + 1),
      q_value_for_plot = pmax(q_value, 1e-50),
      neg_log10_q_capped = q_to_neglog10(q_value, cap = 50),
      q_was_zero = !is.na(q_value) & q_value == 0,
      q_value_raw = ifelse(!is.na(q_value) & q_value == 0, "<1e-300", as.character(q_value)),
      regulation = dplyr::case_when(
        Fold_Change > 0 ~ "Upregula NQO1",
        Fold_Change < 0 ~ "Downregula NQO1",
        TRUE ~ "Sem direção"
      ),
      score = neg_log10_q_capped * abs(Fold_Change)
    ) |>
    dplyr::filter(!is.na(p_value), !is.na(q_value), !is.na(Fold_Change))
}

rank_nqo1_drugs <- function(df) {
  # Unifica nomes de fármacos sem diferenciar maiusculas e minusculas.
  df |>
    dplyr::group_by(Drug_Name_Norm) |>
    dplyr::summarise(
      Drug_Name = names(sort(table(Drug_Name), decreasing = TRUE))[1],
      mean_fc = mean(Fold_Change, na.rm = TRUE),
      median_fc = median(Fold_Change, na.rm = TRUE),
      max_abs_fc = max(abs(Fold_Change), na.rm = TRUE),
      min_p_value = min(p_value, na.rm = TRUE),
      min_q_value = min(q_value, na.rm = TRUE),
      n_tests = dplyr::n(),
      n_sources = dplyr::n_distinct(source_database),
      sources = paste(sort(unique(source_database)), collapse = "; "),
      score = max(score, na.rm = TRUE),
      direction = dplyr::case_when(
        median_fc > 0 ~ "Upregula NQO1",
        median_fc < 0 ~ "Downregula NQO1",
        TRUE ~ "Sem direção"
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(min_q_value, dplyr::desc(score), dplyr::desc(max_abs_fc))
}

plot_volcano <- function(df) {
  top_labels <- df |>
    dplyr::arrange(q_value, dplyr::desc(abs(Fold_Change))) |>
    dplyr::slice_head(n = 12)

  ggplot2::ggplot(df, ggplot2::aes(x = Fold_Change, y = neg_log10_q_capped)) +
    ggplot2::geom_point(ggplot2::aes(color = regulation, shape = q_was_zero), alpha = 0.74, size = 2.1) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    ggplot2::geom_text(
      data = top_labels,
      ggplot2::aes(label = Drug_Name),
      size = 2.8,
      check_overlap = TRUE,
      vjust = -0.65,
      color = "grey20"
    ) +
    ggplot2::scale_color_manual(values = c(
      "Upregula NQO1" = "#1B9E77",
      "Downregula NQO1" = "#D95F02",
      "Sem direção" = "grey60"
    )) +
    ggplot2::scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +
    ggplot2::labs(
      title = "Volcano de fármacos moduladores de NQO1",
      subtitle = "q-values iguais a zero foram capados para evitar escala distorcida",
      x = "Fold change",
      y = "-log10(q-value), capado em 50",
      color = "Regulação",
      shape = "q-value = 0"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_dose_response <- function(df, max_drugs = 9) {
  if (!"Dose" %in% names(df)) df$Dose <- "Dose não informada"
  if (!"Time" %in% names(df)) df$Time <- NA_character_

  dose_df_all <- df |>
    dplyr::mutate(
      Dose = as.character(Dose),
      Time = as.character(Time),
      Dose = dplyr::if_else(is.na(Dose) | Dose == "", "Dose não informada", Dose),
      dose_um = extract_dose_um(Dose),
      time_label = dplyr::if_else(is.na(Time) | Time == "", "Tempo não informado", Time),
      drug_key = Drug_Name_Norm,
      Drug_Display = dplyr::if_else(
        stringr::str_detect(drug_key, "doxo|doxorubic|doxorubicin"),
        "Doxorubicin",
        Drug_Name
      )
    ) |>
    dplyr::filter(!is.na(dose_um), dose_um > 0)

  if (nrow(dose_df_all) == 0) {
    stop("Não ha doses numericas na planilha para montar dose-response.")
  }

  dose_drugs_main <- dose_df_all |>
    dplyr::group_by(Drug_Display) |>
    dplyr::summarise(
      n_doses = dplyr::n_distinct(dose_um),
      n_rows = dplyr::n(),
      min_q_value = min(q_value, na.rm = TRUE),
      score = max(score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(n_doses >= 2) |>
    dplyr::arrange(min_q_value, dplyr::desc(score)) |>
    dplyr::slice_head(n = max_drugs)

  doxo_drug <- dose_df_all |>
    dplyr::filter(Drug_Display == "Doxorubicin") |>
    dplyr::group_by(Drug_Display) |>
    dplyr::summarise(
      n_doses = dplyr::n_distinct(dose_um),
      n_rows = dplyr::n(),
      min_q_value = min(q_value, na.rm = TRUE),
      score = max(score, na.rm = TRUE),
      .groups = "drop"
    )

  dose_drugs <- dplyr::bind_rows(dose_drugs_main, doxo_drug) |>
    dplyr::distinct(Drug_Display, .keep_all = TRUE) |>
    dplyr::arrange(min_q_value, dplyr::desc(score))

  if (nrow(dose_drugs) == 0) {
    dose_drugs <- dose_df_all |>
      dplyr::group_by(Drug_Display) |>
      dplyr::summarise(
        n_rows = dplyr::n(),
        min_q_value = min(q_value, na.rm = TRUE),
        score = max(score, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(min_q_value, dplyr::desc(score)) |>
      dplyr::slice_head(n = max_drugs)
  }

  dose_df <- dose_df_all |>
    dplyr::filter(Drug_Display %in% dose_drugs$Drug_Display) |>
    dplyr::mutate(
      Drug_Display = factor(Drug_Display, levels = rev(dose_drugs$Drug_Display)),
      series = paste(source_database, Cell_Line, time_label, sep = " | ")
    )

  ggplot2::ggplot(dose_df, ggplot2::aes(x = dose_um, y = Fold_Change)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    ggplot2::geom_line(
      ggplot2::aes(group = series, color = source_database),
      alpha = 0.35,
      linewidth = 0.45
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = neg_log10_q_capped, fill = regulation),
      shape = 21,
      color = "grey25",
      stroke = 0.25,
      alpha = 0.88
    ) +
    ggplot2::facet_wrap(~ Drug_Display, scales = "free_y", ncol = 3) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c("CMAP" = "#4C78A8", "LINCS L1000" = "#7F7F7F", "CREEDS" = "#B279A2")) +
    ggplot2::scale_fill_manual(values = c(
      "Upregula NQO1" = "#1B9E77",
      "Downregula NQO1" = "#D95F02",
      "Sem direção" = "grey70"
    )) +
    ggplot2::scale_size_continuous(range = c(1.5, 5.5)) +
    ggplot2::labs(
      title = "Dose-response dos fármacos que modulam NQO1",
      subtitle = "Somente registros com dose numerica; Doxorubicin incluída a partir dos registros LINCS com 10 uM",
      x = "Dose (uM, escala log10)",
      y = "Fold change de NQO1",
      size = "-log10(q) capado",
      color = "Base",
      fill = "Regulação"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

load_pubchem <- function() {
  path <- find_first_existing(c(
    file.path(DIRS$tables_pubchem, "pubchem_associacoes_compostos.csv"),
    file.path(ROOT, "output", "pubchem_associacoes_compostos.csv"),
    file.path(ROOT, "pubchem_associacoes_compostos.csv")
  ), required = FALSE, label = "resultado PubChem")

  if (is.na(path)) {
    audit_log("WARN", "pubchem", "Resultado PubChem não encontrado; gráfico PubChem não será recriado.")
    return(data.frame())
  }
  read_csv_flexible(path)
}

plot_pubchem_counts <- function(pubchem_df) {
  if (nrow(pubchem_df) == 0 || !"Gene" %in% names(pubchem_df)) return(NULL)

  pubchem_df |>
    dplyr::count(Gene, Direction, name = "n_compounds") |>
    dplyr::arrange(dplyr::desc(n_compounds)) |>
    dplyr::slice_head(n = 20) |>
    ggplot2::ggplot(ggplot2::aes(x = n_compounds, y = reorder(Gene, n_compounds), fill = Direction)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::scale_fill_manual(values = c("up" = "#1B9E77", "down" = "#D95F02"), na.value = "grey60") +
    ggplot2::labs(
      title = "Top genes AMPK por compostos associados no PubChem",
      x = "Número de compostos",
      y = NULL,
      fill = "Regulação do gene"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_ppi_network <- function() {
  deg_path <- find_first_existing(c(
    file.path(DIRS$tables_ampk, "expressao_diferencial_via_ampk.csv"),
    file.path(ROOT, "expressao_diferencial_via_ampk.csv")
  ), required = FALSE, label = "expressao_diferencial_via_ampk.csv")
  links_path <- find_first_existing(file.path(DIRS$data_raw, "9606.protein.links.v11.5.min700.txt.gz"),
                                    required = FALSE, label = "STRING links")
  info_path <- find_first_existing(file.path(DIRS$data_raw, "9606.protein.info.v11.5.txt.gz"),
                                   required = FALSE, label = "STRING protein info")

  if (is.na(deg_path) || is.na(links_path) || is.na(info_path)) {
    audit_log("WARN", "ppi", "Arquivos DEG/STRING ausentes; rede PPI não será gerada.")
    return(NULL)
  }

  deg <- read_csv_flexible(deg_path) |>
    dplyr::mutate(gene = as.character(gene))
  genes <- unique(deg$gene)

  audit_log("INFO", "ppi", "Lendo STRING protein info.")
  info <- as.data.frame(data.table::fread(info_path, quote = "", stringsAsFactors = FALSE))
  names(info) <- standardize_names(names(info))
  names(info) <- stringr::str_replace_all(names(info), "^#", "")
  if ("X_string_protein_id" %in% names(info)) {
    names(info)[names(info) == "X_string_protein_id"] <- "string_protein_id"
  }
  info_sub <- info |>
    dplyr::filter(preferred_name %in% genes) |>
    dplyr::select(string_protein_id, preferred_name)

  if (nrow(info_sub) == 0) {
    audit_log("WARN", "ppi", "Nenhum gene AMPK foi mapeado para STRING.")
    return(NULL)
  }

  audit_log("INFO", "ppi", "Lendo STRING links.")
  links <- as.data.frame(data.table::fread(links_path, header = TRUE, stringsAsFactors = FALSE))
  names(links) <- standardize_names(names(links))
  ids <- unique(info_sub$string_protein_id)
  links_sub <- links |>
    dplyr::filter(protein1 %in% ids, protein2 %in% ids) |>
    dplyr::left_join(info_sub, by = c("protein1" = "string_protein_id")) |>
    dplyr::rename(from = preferred_name) |>
    dplyr::left_join(info_sub, by = c("protein2" = "string_protein_id")) |>
    dplyr::rename(to = preferred_name) |>
    dplyr::filter(!is.na(from), !is.na(to), from != to) |>
    dplyr::arrange(dplyr::desc(combined_score))

  if (nrow(links_sub) == 0) {
    audit_log("WARN", "ppi", "Nenhuma interacao STRING entre genes AMPK mapeados.")
    return(NULL)
  }

  utils::write.csv(links_sub, file.path(DIRS$tables_ppi, "ppi_ampk_string_arestas.csv"), row.names = FALSE)

  g <- igraph::graph_from_data_frame(
    d = links_sub |> dplyr::select(from, to, combined_score),
    directed = FALSE
  )
  node_df <- data.frame(gene = igraph::V(g)$name) |>
    dplyr::left_join(deg |> dplyr::select(gene, log2FoldChange, padj, direction), by = "gene")
  igraph::V(g)$degree <- igraph::degree(g)
  igraph::V(g)$log2FoldChange <- node_df$log2FoldChange
  igraph::V(g)$direction <- node_df$direction

  layout <- igraph::layout_with_fr(g, weights = igraph::E(g)$combined_score)
  coords <- data.frame(gene = igraph::V(g)$name, x = layout[, 1], y = layout[, 2]) |>
    dplyr::left_join(node_df, by = "gene") |>
    dplyr::mutate(degree = igraph::degree(g)[gene])

  edge_df <- igraph::as_data_frame(g, what = "edges") |>
    dplyr::left_join(coords |> dplyr::select(gene, x, y), by = c("from" = "gene")) |>
    dplyr::rename(x_from = x, y_from = y) |>
    dplyr::left_join(coords |> dplyr::select(gene, x, y), by = c("to" = "gene")) |>
    dplyr::rename(x_to = x, y_to = y)

  utils::write.csv(coords, file.path(DIRS$tables_ppi, "ppi_ampk_string_nos.csv"), row.names = FALSE)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = x_from, y = y_from, xend = x_to, yend = y_to, linewidth = combined_score),
      color = "grey72",
      alpha = 0.55
    ) +
    ggplot2::geom_point(
      data = coords,
      ggplot2::aes(x = x, y = y, size = degree, fill = log2FoldChange),
      shape = 21,
      color = "grey20",
      stroke = 0.35
    ) +
    ggplot2::geom_text(
      data = coords |> dplyr::arrange(dplyr::desc(degree)) |> dplyr::slice_head(n = 35),
      ggplot2::aes(x = x, y = y, label = gene),
      size = 2.8,
      check_overlap = TRUE,
      vjust = -0.8
    ) +
    ggplot2::scale_fill_gradient2(low = "#D95F02", mid = "white", high = "#1B9E77") +
    ggplot2::scale_linewidth(range = c(0.15, 1.4), guide = "none") +
    ggplot2::labs(
      title = "Rede PPI da via AMPK",
      subtitle = "STRING Homo sapiens, interacoes min700; cor = log2FC em DCM vs controle",
      x = NULL,
      y = NULL,
      fill = "log2FC",
      size = "Grau"
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  save_plot(p, file.path(DIRS$figures_ppi, "rede_ppi_ampk.png"), width = 10, height = 8)
  audit_log("INFO", "ppi", paste("Rede PPI:", igraph::vcount(g), "nos e", igraph::ecount(g), "arestas."))
  invisible(list(graph = g, nodes = coords, edges = links_sub))
}

cleanup_unneeded_outputs <- function() {
  unneeded <- c(
    file.path(DIRS$figures_nqo1, "boxplot_top_drugs_NQO1.png"),
    file.path(DIRS$figures_nqo1, "TOP20_FARMACOS_NQO1.png"),
    file.path(DIRS$figures_ampk, "PCA_final_publication.png"),
    file.path(DIRS$figures_ampk, "PCA_microarray.png"),
    file.path(DIRS$figures_ampk, "PCA_scree_plot.png")
  )
  for (path in unneeded[file.exists(unneeded)]) {
    file.remove(path)
    audit_log("INFO", "cleanup", paste("Removido por não ser mais necessario:", path))
  }
}

run_downstream_analysis <- function() {
  audit_log("INFO", "downstream", "Início da análise downstream (PASSO 3).")
  cleanup_unneeded_outputs()

  xlsx_path <- find_first_existing(c(
    file.path(DIRS$tables_nqo1, "dgb_resultados_nqo1.xlsx"),
    file.path(ROOT, "dgb_resultados_nqo1.xlsx")
  ), label = "dgb_resultados_nqo1.xlsx")

  dgb <- read_nqo1_xlsx(xlsx_path)
  audit_log("INFO", "nqo1", paste("Linhas validas na planilha NQO1:", nrow(dgb)))
  audit_log("INFO", "nqo1", paste("Fármacos únicos (pos-unificacao case-insensitive):",
                                   dplyr::n_distinct(dgb$Drug_Name_Norm)))

  ranking <- rank_nqo1_drugs(dgb)
  top20 <- ranking |> dplyr::slice_head(n = 20)

  utils::write.csv(dgb, file.path(DIRS$tables_nqo1, "dgb_resultados_nqo1_consolidados.csv"), row.names = FALSE)
  utils::write.csv(ranking, file.path(DIRS$tables_nqo1, "ranking_de_farmacos_nqo1.csv"), row.names = FALSE)
  utils::write.csv(top20, file.path(DIRS$tables_nqo1, "top20_farmacos_nqo1.csv"), row.names = FALSE)

  # Tabela 1: distribuição por plataforma e direção.
  tab1_platform <- ranking |>
    dplyr::mutate(
      platform_group = dplyr::case_when(
        sources == "LINCS L1000" ~ "LINCS L1000",
        sources == "CREEDS" ~ "CREEDS",
        sources == "CMAP" ~ "CMAP",
        TRUE ~ "Multi-plataforma"
      )
    ) |>
    dplyr::count(platform_group, direction, name = "n_farmacos") |>
    tidyr::pivot_wider(
      names_from = direction,
      values_from = n_farmacos,
      values_fill = 0
    )
  utils::write.csv(tab1_platform, file.path(DIRS$tables_nqo1, "farmacos_por_plataforma_e_direcao.csv"), row.names = FALSE)
  audit_log("INFO", "table1", "Tabela 1 gerada: plataformas x direção")

  doxo <- ranking |>
    dplyr::filter(stringr::str_detect(Drug_Name_Norm, "doxo|doxorubic|doxorubicin"))
  utils::write.csv(doxo, file.path(DIRS$tables_nqo1, "resumo_doxorrubicina_nqo1.csv"), row.names = FALSE)
  if (nrow(doxo) > 0) {
    audit_log("INFO", "nqo1", paste("Doxorubicin encontrada (pos-unificacao); melhor q-value:",
                                     signif(min(doxo$min_q_value), 4)))
    audit_log("INFO", "nqo1", paste("Direção da doxorubicina:", doxo$direction[1],
                                     "| mean_fc:", signif(doxo$mean_fc[1], 4)))
  } else {
    audit_log("WARN", "nqo1", "Doxorubicin não foi encontrada no ranking.")
  }

  save_plot(plot_volcano(dgb), file.path(DIRS$figures_nqo1, "volcano_farmacos_nqo1.png"), 10, 7.5)
  save_plot(plot_dose_response(dgb), file.path(DIRS$figures_nqo1, "curva_dose_resposta_farmacos_nqo1.png"), 11.5, 8.5)

  pubchem <- load_pubchem()
  if (nrow(pubchem) > 0) {
    utils::write.csv(pubchem, file.path(DIRS$tables_pubchem, "pubchem_associacoes_compostos.csv"), row.names = FALSE)
    pubchem_named <- pubchem |> dplyr::filter(!is.na(Name), Name != "")
    utils::write.csv(pubchem_named, file.path(DIRS$tables_pubchem, "pubchem_associacoes_compostos_com_nomes.csv"), row.names = FALSE)
    p_pub <- plot_pubchem_counts(pubchem)
    if (!is.null(p_pub)) {
      save_plot(p_pub, file.path(DIRS$figures_pubchem, "pubchem_compostos_por_gene.png"), 9, 6.5)
    }
  }

  plot_ppi_network()

  n_down <- sum(ranking$direction == "Downregula NQO1")
  n_up <- sum(ranking$direction == "Upregula NQO1")
  n_multi <- sum(stringr::str_detect(ranking$sources, ";"))
  n_lincs <- sum(stringr::str_detect(ranking$sources, "LINCS"))
  n_creeds <- sum(stringr::str_detect(ranking$sources, "CREEDS"))
  n_cmap <- sum(stringr::str_detect(ranking$sources, "CMAP"))

  report <- c(
    "# Auditoria NQO1/PubChem/PPI",
    "",
    paste("- Data:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("- Raiz:", ROOT),
    paste("- Fonte única dos fármacos:", xlsx_path),
    paste("- Abas processadas:", paste(openxlsx::getSheetNames(xlsx_path), collapse = ", ")),
    paste("- Linhas NQO1 validas:", nrow(dgb)),
    paste("- Fármacos únicos (pos-unificacao case-insensitive):", nrow(ranking)),
    paste("- Downreguladores:", n_down, "(", round(100 * n_down / nrow(ranking), 1), "%)"),
    paste("- Upreguladores:", n_up, "(", round(100 * n_up / nrow(ranking), 1), "%)"),
    paste("- LINCS L1000:", n_lincs, "(", round(100 * n_lincs / nrow(ranking), 1), "%)"),
    paste("- CREEDS:", n_creeds, "(", round(100 * n_creeds / nrow(ranking), 1), "%)"),
    paste("- CMAP:", n_cmap, "(", round(100 * n_cmap / nrow(ranking), 1), "%)"),
    paste("- Multi-plataforma:", n_multi, "(", round(100 * n_multi / nrow(ranking), 1), "%)"),
    paste("- Top 1:", top20$Drug_Name[1]),
    paste("- Doxorubicin no ranking:", ifelse(nrow(doxo) > 0, "sim", "não")),
    if (nrow(doxo) > 0) {
      paste("- Doxorubicin direção:", doxo$direction[1], "| mean_fc:", signif(doxo$mean_fc[1], 4),
            "| min_q:", signif(doxo$min_q_value[1], 4))
    } else {
      "- Doxorubicin: não encontrada"
    },
    "- Volcano revisado: q-values iguais a zero capados em -log10(q) = 50.",
    "- Boxplot top drugs removido.",
    "- PCA removido.",
    paste("- Tabela top20 preservada como tabela:", file.path(DIRS$tables_nqo1, "top20_farmacos_nqo1.csv")),
    "- Gráfico TOP20_FARMACOS_NQO1 removido.",
    "- Unificacao case-insensitive aplicada no ranking.",
    paste("- Figuras NQO1:", DIRS$figures_nqo1),
    paste("- Figura PPI:", file.path(DIRS$figures_ppi, "rede_ppi_ampk.png")),
    paste("- Log estruturado:", file.path(DIRS$audit, "eventos_de_auditoria.csv")),
    "",
    "## Top 20",
    paste(capture.output(print(top20[, c("Drug_Name", "min_q_value", "mean_fc", "score", "direction", "sources")],
                               row.names = FALSE)), collapse = "\n")
  )
  writeLines(report, file.path(DIRS$audit, "relatorio_nqo1_pubchem_ppi.md"), useBytes = TRUE)

  audit_log("INFO", "downstream", "Análise downstream concluída.")
  invisible(NULL)
}

# ==============================================================================
# PARTE C. Auditoria de dados
# ==============================================================================
# Checagens de integridade, valores ausentes, duplicatas, outliers e
# consistência entre amostras e metadados. O resultado e gravado em
# output/audit/auditoria_de_dados_resumo.txt.
# ==============================================================================

audit_pct_na <- function(df) {
  vapply(df, function(x) round(100 * mean(is.na(x)), 2), numeric(1))
}

audit_z_outliers <- function(x, thr = 3.0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(unique(x)) <= 1) return(0L)
  mu <- mean(x)
  sd_pop <- sqrt(mean((x - mu)^2))
  if (sd_pop == 0) return(0L)
  as.integer(sum(abs((x - mu) / sd_pop) > thr))
}

audit_csv_stats <- function(path) {
  if (!file.exists(path)) return("ausente")
  df <- tryCatch(
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) return("erro de leitura")
  paste0(nrow(df), " linhas x ", ncol(df), " colunas")
}

run_data_audit <- function() {
  audit_log("INFO", "audit", "Início da auditoria de dados.")
  lines <- c(
    "# Auditoria de dados",
    "",
    paste("Data:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    ""
  )

  # 1. Tabela NQO1 limpa.
  dgb_path <- file.path(DIRS$tables_nqo1, "dgb_resultados_nqo1_consolidados.csv")
  if (file.exists(dgb_path)) {
    dgb <- read.csv(dgb_path, check.names = FALSE, stringsAsFactors = FALSE)
    na_pct <- audit_pct_na(dgb)
    lines <- c(
      lines,
      "## 1. Tabela NQO1 limpa",
      paste("Dimensoes:", nrow(dgb), "linhas x", ncol(dgb), "colunas"),
      paste("Linhas duplicadas:", sum(duplicated(dgb))),
      paste("Coluna-indice X1 presente:", "X1" %in% names(dgb)),
      paste("p_value igual a zero:", sum(dgb$p_value == 0, na.rm = TRUE)),
      paste("q_value igual a zero:", sum(dgb$q_value == 0, na.rm = TRUE)),
      paste("Colunas com mais de 5% de valores ausentes:",
            if (any(na_pct > 5)) paste(names(na_pct)[na_pct > 5], collapse = ", ") else "nenhuma"),
      ""
    )

    # 2. Doxorrubicina.
    norm <- ifelse(is.na(dgb$Drug_Name_Norm), "", dgb$Drug_Name_Norm)
    doxo <- dgb[grepl("doxo", norm), , drop = FALSE]
    if (nrow(doxo) > 0) {
      lines <- c(
        lines,
        "## 2. Doxorubicina e NQO1",
        paste("Registros:", nrow(doxo)),
        paste("Fold Change medio:", round(mean(doxo$Fold_Change, na.rm = TRUE), 4)),
        paste("Fold Change mediano:", round(median(doxo$Fold_Change, na.rm = TRUE), 4)),
        paste("Maior |Fold Change|:", round(max(abs(doxo$Fold_Change), na.rm = TRUE), 4)),
        paste("Menor q-value:", min(doxo$q_value, na.rm = TRUE)),
        ""
      )
    }
  }

  # 3. Expressão diferencial completa e NQO1.
  deg_path <- file.path(DIRS$tables_ampk, "expressao_diferencial_completa.csv")
  if (file.exists(deg_path)) {
    deg <- read.csv(deg_path, check.names = FALSE, stringsAsFactors = FALSE)
    nqo1 <- deg[deg$gene == "NQO1", , drop = FALSE]
    lines <- c(
      lines,
      "## 3. Expressão diferencial",
      paste("Genes:", nrow(deg), "| genes duplicados:", sum(duplicated(deg$gene))),
      paste("Genes significativos (padj < 0.05):", sum(deg$significant, na.rm = TRUE))
    )
    if (nrow(nqo1) > 0) {
      lines <- c(
        lines,
        paste("NQO1 log2FC:", round(nqo1$log2FoldChange[1], 4),
              "| padj:", signif(nqo1$padj[1], 4),
              "| direção:", nqo1$direction[1]),
        paste("NQO1 supera o limiar |log2FC| > 0.585:", abs(nqo1$log2FoldChange[1]) > 0.585)
      )
    }
    lines <- c(lines, "")
  }

  # 4. Outliers por Z-score.
  if (file.exists(dgb_path)) {
    cols <- c("p_value", "q_value", "Fold_Change", "Specificity", "score", "neg_log10_q_capped")
    cols <- cols[cols %in% names(dgb)]
    z <- vapply(cols, function(col) audit_z_outliers(dgb[[col]]), integer(1))
    lines <- c(
      lines,
      "## 4. Outliers (|Z| > 3)",
      paste(names(z), z, sep = ": ", collapse = " | "),
      ""
    )
  }

  # 5. Matriz RPKM e NQO1.
  rpkm_path <- file.path(DIRS$data_raw, "GSE116250_rpkm.txt.gz")
  if (file.exists(rpkm_path)) {
    con <- gzfile(rpkm_path, "rt")
    header <- strsplit(readLines(con, n = 1, warn = FALSE), "\t")[[1]]
    samples <- header[-c(1, 2)]
    n_missing <- 0L
    total_cells <- 0L
    n_genes <- 0L
    nqo1_vals <- numeric(0)
    repeat {
      line <- readLines(con, n = 1, warn = FALSE)
      if (length(line) == 0) break
      parts <- strsplit(line, "\t")[[1]]
      if (length(parts) < 3) next
      n_genes <- n_genes + 1L
      vals <- parts[-c(1, 2)]
      n_missing <- n_missing + sum(trimws(vals) == "")
      total_cells <- total_cells + length(vals)
      if (parts[2] == "NQO1") nqo1_vals <- suppressWarnings(as.numeric(vals[trimws(vals) != ""]))
    }
    close(con)
    classes <- sub("^([A-Za-z]+).*$", "\\1", samples)
    class_counts <- table(classes)
    lines <- c(
      lines,
      "## 5. Matriz RPKM",
      paste("Amostras:", length(samples), "| genes:", n_genes),
      paste("Classes:", paste(names(class_counts), as.integer(class_counts), sep = " = ", collapse = ", ")),
      paste("Celulas ausentes:", n_missing, "de", total_cells,
            "(", round(100 * n_missing / total_cells, 4), "%)"),
      paste("NQO1 presente:", length(nqo1_vals) > 0)
    )
    if (length(nqo1_vals) > 0) {
      lines <- c(
        lines,
        paste("NQO1 RPKM medio:", round(mean(nqo1_vals), 3),
              "| mediano:", round(median(nqo1_vals), 3))
      )
    }
    lines <- c(lines, "")
  }

  # 6. Consistência entre amostras e metadados.
  pheno_path <- file.path(DIRS$tables_ampk, "metadados_das_amostras.csv")
  if (file.exists(rpkm_path) && file.exists(pheno_path)) {
    pheno <- read.csv(pheno_path, check.names = FALSE, stringsAsFactors = FALSE)
    rpk_set <- unique(as.character(samples))
    pheno_set <- unique(as.character(pheno$title))
    lines <- c(
      lines,
      "## 6. Consistencia amostra x metadado",
      paste("Amostras na matriz:", length(rpk_set), "| amostras no metadado:", length(pheno_set)),
      paste("Somente na matriz:", paste(setdiff(rpk_set, pheno_set), collapse = ", ")),
      paste("Somente no metadado:", paste(setdiff(pheno_set, rpk_set), collapse = ", ")),
      ""
    )
  }

  # 7. Duplicatas na rede PPI.
  edges_path <- file.path(DIRS$tables_ppi, "ppi_ampk_string_arestas.csv")
  nodes_path <- file.path(DIRS$tables_ppi, "ppi_ampk_string_nos.csv")
  if (file.exists(edges_path) && file.exists(nodes_path)) {
    edges <- read.csv(edges_path, check.names = FALSE, stringsAsFactors = FALSE)
    nodes <- read.csv(nodes_path, check.names = FALSE, stringsAsFactors = FALSE)
    lines <- c(
      lines,
      "## 7. Rede PPI",
      paste("Arestas:", nrow(edges), "| duplicadas:", sum(duplicated(edges))),
      paste("Nos:", nrow(nodes), "| duplicados:", sum(duplicated(nodes))),
      ""
    )
  }

  # 8. GSEA da via AMPK.
  gsea_path <- file.path(DIRS$tables_ampk, "gsea_via_ampk.csv")
  if (file.exists(gsea_path)) {
    gsea <- read.csv(gsea_path, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(gsea) > 0) {
      lines <- c(
        lines,
        "## 8. GSEA da via AMPK",
        paste("NES:", round(gsea$NES[1], 4), "| padj:", round(gsea$padj[1], 4),
              "| significativo (padj < 0.05):", gsea$padj[1] < 0.05),
        ""
      )
    }
  }

  # 9. Integridade dos outputs finais.
  outputs <- c(
    expressao_diferencial_completa.csv = file.path(DIRS$tables_ampk, "expressao_diferencial_completa.csv"),
    expressao_nqo1_diferencial.csv = file.path(DIRS$tables_ampk, "expressao_nqo1_diferencial.csv"),
    expressao_diferencial_via_ampk.csv = file.path(DIRS$tables_ampk, "expressao_diferencial_via_ampk.csv"),
    gsea_via_ampk.csv = file.path(DIRS$tables_ampk, "gsea_via_ampk.csv"),
    gsva_via_ampk.csv = file.path(DIRS$tables_ampk, "gsva_via_ampk.csv"),
    resumo_doxorrubicina_nqo1.csv = file.path(DIRS$tables_nqo1, "resumo_doxorrubicina_nqo1.csv"),
    dgb_resultados_nqo1.xlsx = file.path(DIRS$tables_nqo1, "dgb_resultados_nqo1.xlsx")
  )
  lines <- c(lines, "## 9. Integridade dos outputs finais")
  for (name in names(outputs)) {
    path <- outputs[[name]]
    if (grepl("\\.xlsx$", name)) {
      lines <- c(lines, paste(name, ":", if (file.exists(path)) "presente" else "ausente"))
    } else {
      lines <- c(lines, paste(name, ":", audit_csv_stats(path)))
    }
  }

  writeLines(lines, file.path(DIRS$audit, "auditoria_de_dados_resumo.txt"), useBytes = TRUE)
  audit_log("INFO", "audit", "Auditoria concluída; resumo em output/audit/auditoria_de_dados_resumo.txt")
  invisible(NULL)
}

# ==============================================================================
# Execução principal
# ==============================================================================

main <- function() {
  audit_log("INFO", "start", paste("Raiz do projeto:", ROOT))
  run_upstream_analysis()
  run_downstream_analysis()
  run_data_audit()

  utils::write.csv(AUDIT_EVENTS, file.path(DIRS$audit, "eventos_de_auditoria.csv"), row.names = FALSE)
  audit_log("INFO", "done", "Análise concluída.")
}

main()
