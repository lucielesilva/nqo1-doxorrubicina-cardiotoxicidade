#!/usr/bin/env Rscript

# ==============================================================================
# 00 — CONFIGURAÇÃO E DADOS (a "base" de todos os scripts)
# ==============================================================================
#
# PARA QUE SERVE ESTE ARQUIVO?
#   Ele NÃO faz análise. Ele apenas:
#     1) descobre em que pasta o projeto está;
#     2) cria as pastas de dados e de resultados;
#     3) baixa os dados do GEO (uma vez) e guarda um "cache" para não baixar
#        tudo de novo a cada aula;
#     4) baixa a lista de genes da via AMPK do KEGG (também com cache);
#     5) define funções e caminhos que os outros scripts vão usar.
#
# POR QUE SEPARAR ISSO?
#   Os scripts 01 a 06 precisam dos MESMOS dados. Se cada um baixasse
#   tudo de novo, a aula ficaria lenta e repetitiva. Aqui a gente baixa uma vez
#   e reutiliza. É como preparar os ingredientes antes de cozinhar.
#
# COMO USAR?
#   Normalmente você NÃO precisa rodar este arquivo sozinho: os scripts 01–06
#   carregam ele automaticamente com source(). Mas, se quiser baixar tudo antes
#   da aula, rode (a partir da raiz do projeto):
#
#     & 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' 'tutoriais\00_configuracao_e_dados.R'
#
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(42)   # deixa resultados aleatórios reprodutíveis (layout, permutações)

# ------------------------------------------------------------------------------
# 0.1) ONDE ESTE SCRIPT ESTÁ? (para funcionar de qualquer pasta)
# ------------------------------------------------------------------------------
# commandArgs() devolve os argumentos usados para chamar o R. Quando rodamos
# com Rscript, existe um argumento "--file=caminho/do/script.R". Usamos isso
# para descobrir a pasta do script. Se o arquivo for carregado com source() a
# partir do RStudio, usamos a pasta da aba aberta. Se nada funcionar, usamos a
# pasta atual.
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

# Os scripts ficam em tutoriais/, mas os dados e os resultados ficam na RAIZ do
# projeto (data/ e output/). Por isso subimos um nível: procuramos, a partir da
# pasta do script, a primeira pasta que tenha o DESCRIPTION do projeto (ou as
# pastas data/ e output/). Assim funciona tanto rodando da raiz quanto de dentro
# de tutoriais/, sem espalhar cópias de data/ e output/ por aí.
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
  pasta_script   # último recurso: assume que o script já está na raiz
}

ROOT <- detectar_raiz(pasta_do_script())  # raiz do projeto (a pasta ACIMA de tutoriais/)
setwd(ROOT)                 # todos os caminhos relativos partem daqui

# ------------------------------------------------------------------------------
# 0.2) PASTAS DO PROJETO
# ------------------------------------------------------------------------------
# data/raw       -> arquivos brutos baixados do GEO e cache do STRING
# data/processed -> objetos .rds já tratados (cache de dados e do limma)
# output/tables  -> tabelas .csv para abrir no Excel
# output/figures -> figuras .png para slides
RAW          <- file.path(ROOT, "data", "raw")
PROCESSED    <- file.path(ROOT, "data", "processed")
STRING_CACHE <- file.path(RAW, "string_cache")   # cache de downloads do STRING
TABLES       <- file.path(ROOT, "output", "tables")
FIGURES      <- file.path(ROOT, "output", "figures")

# dir.create() cria a pasta só se ela ainda não existir.
for (d in c(RAW, PROCESSED, STRING_CACHE, TABLES, FIGURES)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------------------------
# 0.3) IDs DOS DADOS
# ------------------------------------------------------------------------------
GSE_ID  <- "GSE116250"   # estudo público do GEO: RNA-seq de ventrículo esquerdo
AMPK_ID <- "hsa04152"    # via KEGG "AMPK signaling pathway" em humano (hsa)

# ------------------------------------------------------------------------------
# 0.4) PACOTES NECESSÁRIOS PARA BAIXAR/ORGANIZAR OS DADOS
# ------------------------------------------------------------------------------
# Cada script de análise carrega também os pacotes específicos dele (limma,
# fgsea, GSVA etc.). Aqui carregamos só o básico para os dados.
suppressPackageStartupMessages({
  library(GEOquery)   # baixar dados do GEO
  library(Biobase)    # acessar metadados (pData)
  library(dplyr)      # manipular tabelas
  library(readr)      # ler e escrever CSV
})

# ------------------------------------------------------------------------------
# 0.5) FUNÇÕES AUXILIARES
# ------------------------------------------------------------------------------
# Mensagem de "passo" no console. Ajuda a acompanhar o que está acontecendo.
step <- function(...) message("== ", ...)

# Salva um gráfico ggplot com padrão de publicação (300 dpi, fundo branco).
# p = objeto do gráfico; file = nome do arquivo; w/h = largura/altura em polegadas.
save_plot <- function(p, file, w = 9, h = 6) {
  ggplot2::ggsave(file.path(FIGURES, file), p,
                  width = w, height = h, dpi = 300, bg = "white")
}

# ==============================================================================
# 0.6) FUNÇÃO: preparar_dados()
# ==============================================================================
# O QUE ELA FAZ?
#   1) Se já existe o cache "dados_GSE116250.rds", carrega e devolve.
#   2) Senão, baixa do GEO:
#        - metadados das amostras (doença, sexo, idade);
#        - matriz de expressão em RPKM (gene x amostra).
#   3) Organiza tudo, salva o cache e devolve uma lista com:
#        expr  = matriz de expressão (linhas = genes, colunas = amostras)
#        pheno = tabela de informações das amostras
#
# GLOSSÁRIO RÁPIDO:
#   RPKM = Reads Per Kilobase per Million mapped reads. É uma medida de
#          expressão já normalizada pelo tamanho do gene e pela profundidade
#          do sequenciamento. Quanto maior, mais aquele gene foi "lido".
#   GEQ   = repositório público de dados de expressão (NCBI GEO).
# ==============================================================================
preparar_dados <- function(forcar = FALSE) {
  cache <- file.path(PROCESSED, paste0("dados_", GSE_ID, ".rds"))

  # Se o cache existe e o usuário não pediu para forçar, usa o cache.
  if (file.exists(cache) && !forcar) {
    message("Dados já preparados. Carregando cache: ", cache)
    return(readRDS(cache))
  }

  # ---------------------------------------------------------------------------
  # PASSO A: baixar os METADADOS do estudo
  # ---------------------------------------------------------------------------
  step("Baixando metadados de ", GSE_ID, " no GEO")
  eset <- GEOquery::getGEO(GSE_ID, GSEMatrix = TRUE, getGPL = FALSE,
                           destdir = RAW)[[1]]
  meta <- Biobase::pData(eset)   # tabela com uma linha por amostra

  # Montamos só as colunas que interessam. Os nomes "Sex:ch1", "age:ch1" etc.
  # são os nomes que o GEO usa para esses campos neste estudo.
  pheno <- data.frame(
    sample = as.character(meta$title),                    # ex.: NF1, DCM8...
    sex    = as.character(meta[["Sex:ch1"]]),             # covariável
    age    = suppressWarnings(as.numeric(meta[["age:ch1"]])), # covariável
    stringsAsFactors = FALSE
  )

  # Traduzimos os rótulos originais para três grupos curtos:
  #   CTRL = coração que não falhou (non-failing)
  #   DCM  = cardiomiopatia dilatada
  #   ICM  = cardiomiopatia isquêmica (outra doença; usada só como contexto no UMAP)
  pheno$group <- dplyr::case_when(
    meta[["disease:ch1"]] == "non-failing"             ~ "CTRL",
    meta[["disease:ch1"]] == "dilated cardiomyopathy"  ~ "DCM",
    meta[["disease:ch1"]] == "ischemic cardiomyopathy" ~ "ICM",
    TRUE ~ NA_character_
  )

  # Remove amostras sem grupo identificado.
  pheno <- pheno[!is.na(pheno$group), , drop = FALSE]

  # ---------------------------------------------------------------------------
  # PASSO B: baixar a MATRIZ DE EXPRESSÃO em RPKM
  # ---------------------------------------------------------------------------
  # O arquivo suplementar do GEO tem a expressão de cada gene em cada amostra.
  # Se já baixamos antes, não baixamos de novo.
  rpkm_path <- file.path(RAW, paste0(GSE_ID, "_rpkm.txt.gz"))

  if (!file.exists(rpkm_path)) {
    step("Baixando a matriz de expressão RPKM (pode demorar alguns minutos)")
    GEOquery::getGEOSuppFiles(GSE_ID, baseDir = RAW, fetch_files = TRUE)

    candidatos <- list.files(file.path(RAW, GSE_ID),
                             pattern = "rpkm.*(txt|gz)$", full.names = TRUE)
    if (length(candidatos) == 0) {
      stop("Não encontrei o arquivo RPKM dentro de data/raw/", GSE_ID,
           ". Verifique a internet ou baixe manualmente do GEO.", call. = FALSE)
    }
    file.copy(candidatos[1], rpkm_path, overwrite = TRUE)
  }

  # Lemos o arquivo. As duas primeiras colunas são anotação (Gene e
  # Common_name); o resto são as amostras. Forçamos as colunas de anotação a
  # serem texto para evitar que símbolos como "1-Mar" virem data.
  rpkm <- as.data.frame(readr::read_tsv(
    rpkm_path,
    show_col_types = FALSE,
    col_types = readr::cols(
      .default     = readr::col_double(),
      Gene         = readr::col_character(),
      Common_name  = readr::col_character()
    )
  ))

  # Matriz: linhas = genes (símbolo Common_name), colunas = amostras.
  expr <- as.matrix(rpkm[, -(1:2), drop = FALSE])
  rownames(expr) <- make.unique(rpkm$Common_name)  # duplicados viram ".1"
  storage.mode(expr) <- "double"                   # garante que é número

  # ---------------------------------------------------------------------------
  # PASSO C: alinhar metadados e matriz
  # ---------------------------------------------------------------------------
  # Só ficamos com amostras que aparecem nos DOIS arquivos.
  keep <- intersect(pheno$sample, colnames(expr))
  if (length(keep) == 0) {
    stop("Nenhuma amostra em comum entre metadados e matriz de expressão.",
         call. = FALSE)
  }

  expr  <- expr[, keep, drop = FALSE]
  pheno <- pheno[match(keep, pheno$sample), , drop = FALSE]
  rownames(pheno) <- NULL

  # Salvamos a tabela de metadados para consulta no Excel.
  readr::write_csv(pheno, file.path(TABLES, "metadados_das_amostras.csv"))

  message("Matriz: ", nrow(expr), " genes x ", ncol(expr), " amostras | ",
          "CTRL=", sum(pheno$group == "CTRL"), " ",
          "DCM=",  sum(pheno$group == "DCM"),  " ",
          "ICM=",  sum(pheno$group == "ICM"))

  # ---------------------------------------------------------------------------
  # PASSO D: salvar o cache e devolver
  # ---------------------------------------------------------------------------
  dados <- list(expr = expr, pheno = pheno, gse = GSE_ID)
  saveRDS(dados, cache)
  message("Cache salvo em: ", cache)

  dados
}

# ==============================================================================
# 0.7) FUNÇÃO: baixar_genes_ampk()
# ==============================================================================
# O QUE ELA FAZ?
#   Baixa do KEGG a lista de genes da via "AMPK signaling pathway" (hsa04152)
#   e devolve uma tabela com:
#     entrez_id   = ID interno do KEGG
#     symbol      = símbolo do gene (ex.: PRKAA1)
#     description = descrição
#
# POR QUE AMPK?
#   É uma via ligada a estresse energético, autofagia e sobrevivência celular
#   no coração. Queremos testar se a via COMO UM TODO está perturbada na DCM.
#
# ATENÇÃO: o campo GENE do KEGG vem "plano": posições ímpares = ID Entrez,
# posições pares = "SÍMBOLO; descrição". Por isso separamos de 2 em 2.
# ==============================================================================
baixar_genes_ampk <- function(forcar = FALSE) {
  ampk_path <- file.path(TABLES, "genes_da_via_ampk_kegg.csv")

  # Cache: se já baixamos a via, não precisamos de internet de novo.
  if (file.exists(ampk_path) && !forcar) {
    message("Via AMPK já baixada. Carregando cache: ", ampk_path)
    return(readr::read_csv(ampk_path, show_col_types = FALSE))
  }

  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop("O pacote KEGGREST não está instalado. Instale com:\n",
         "  BiocManager::install('KEGGREST')", call. = FALSE)
  }

  step("Baixando a via AMPK (", AMPK_ID, ") do KEGG")
  kegg <- tryCatch(KEGGREST::keggGet(AMPK_ID)[[1]]$GENE,
                   error = function(e) NULL)

  if (is.null(kegg)) {
    stop("Não foi possível baixar a via AMPK do KEGG. ",
         "Verifique a conexão com a internet e tente novamente.", call. = FALSE)
  }

  defs <- kegg[seq(2, length(kegg), 2)]   # posições pares: "SÍMBOLO; descrição"

  ampk <- data.frame(
    entrez_id   = kegg[seq(1, length(kegg), 2)],      # posições ímpares: Entrez
    symbol      = trimws(sub(";.*$", "", defs)),      # texto antes do ";"
    description = trimws(defs),
    stringsAsFactors = FALSE
  )

  ampk <- ampk[ampk$symbol != "", , drop = FALSE]      # remove símbolos vazios
  ampk <- ampk[!duplicated(ampk$symbol), , drop = FALSE] # remove duplicados

  readr::write_csv(ampk, ampk_path)
  message("Via AMPK: ", nrow(ampk), " genes")
  ampk
}

# ==============================================================================
# 0.8) SE ESTE ARQUIVO FOR RODADO SOZINHO, JÁ PREPARA TUDO
# ==============================================================================
# sys.nframe() == 0 significa "este é o script principal, não fui carregado
# por outro script com source()". Nesse caso, baixamos os dados e a via AMPK.
# Quando o 00 é carregado pelo 01/02/03/04, isso NÃO roda automaticamente.
if (sys.nframe() == 0) {
  dados <- preparar_dados()
  ampk  <- baixar_genes_ampk()
  step("Base pronta! Agora escolha o script da aula: 01_UMAP.R, 02_volcano.R, etc.")
}
