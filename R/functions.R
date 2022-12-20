#' count_matrix_loader
#'
#' @param file_type your raw data file type (presently optimized for txt)
#'
#' @return a count matrix containing raw count values for the experiment with group annotations.

count_matrix_assembly <- function(file_type) {
  count_file <- fs::dir_ls(here::here("data-raw/"),
    regexp = file_type,
    recurse = TRUE
  )
  count_matrix <- openxlsx::read.xlsx(xlsxFile = count_file, sheet = 1)
  count_matrix$Geneid <- stringr::str_remove_all(count_matrix$Geneid, "\\..*")
  rownames(count_matrix) <- count_matrix$Geneid
  count_matrix <- count_matrix %>%
    dplyr::select(-Geneid)

  return(count_matrix)
}

#' Load metadata and sort them according to count matrix input
#'
#' @param file_name the name of teh metadata file (default "metadata.csv")
#'
#' @return metadata file sorted according to count matrix order

load_metadata <- function(file_name) {
  data_file <- fs::dir_ls(here::here("data-raw/"),
    regexp = file_name,
    recurse = T
  )
  metadata <- openxlsx::read.xlsx(xlsxFile = data_file)
  metadata <- metadata |>
    dplyr::mutate(Group = dplyr::case_when(
      Group == "WT_L" ~ "L",
      Group == "WT_CS" ~ "CS",
      Group == "WT_PH" ~ "PH"
    )) |>
    dplyr::filter(Genotype == "WT")
  return(metadata)
}

#' Remove identified outlier samples from metadata
#'
#' @param outlier
#'
#' @return

removeOutlierMeta <- function(metadata, outlier) {
  metadata <- metadata |>
    dplyr::filter(!Sample == outlier)
}

#' Remove identified outliers from counts and arrange in same order of metadata
#'
#' @param counts
#' @param metadata
#' @param outlier
#'
#' @return

removeOutlierCounts <- function(counts, metadata, outlier) {
  counts <- counts %>%
    dplyr::select(-outlier) |>
    dplyr::select(metadata$Sample)
}

#' Generate design matrix
#'
#' @param metadata a metadata object generated through the load_metadata function
#'
#' @return a design matrix file


Generate_design_matrix <- function(metadata) {
  design <- stats::model.matrix(~ 0 + Group, metadata)
  colnames(design) <-
    stringr::str_remove_all(colnames(design), "\\(|\\)|Group|:")
  return(design)
}

#' Generate contrast matrix
#'
#' @param design
#'
#' @return

Generate_ctrst <- function(design) {
  limma::makeContrasts(
    Liv_vs_CS = L - CS,
    Liv_vs_PH = L - PH,
    CS_vs_PH = CS - PH,
    levels = design
  )
}

#'  DGE analysis of counts using generated contrast matrix
#'
#' @param counts
#' @param metadata
#' @param design
#' @param ctrsts
#'
#' @return

DGEanalysis <- function(counts, metadata, design, ctrsts) {
  group <- as.matrix(metadata$Group)
  RNAseq <- edgeR::DGEList(counts = counts,
                           group = group)
  keep <- edgeR::filterByExpr(RNAseq,
                              design = design)
  RNAseq <- RNAseq[keep, , keep.lib.sizes = F]
  RNAseq <- edgeR::calcNormFactors(RNAseq)
  cpm_matrix <- edgeR::cpm(RNAseq,
                           log = T)
  key <- clusterProfiler::bitr(rownames(cpm_matrix),
                               fromType = "ENSEMBL",
                               toType = "SYMBOL",
                               OrgDb = "org.Mm.eg.db")
  RNAseq <- edgeR::estimateDisp(RNAseq,
                                design)
  efit <- edgeR::glmQLFit(RNAseq,
                          design)
  dgeResults <- apply(ctrsts, 2, . %>%
    edgeR::glmQLFTest(glmfit = efit, contrast = .) %>%
    edgeR::topTags(n = Inf, p.value = 1) %>%
    magrittr::extract2("table") %>%
    data.table::as.data.table(keep.rownames = TRUE))
  dgeResults_annotated <- dgeResults
  for (i in 1:length(dgeResults_annotated)) {
    data.table::setnames(dgeResults_annotated[[i]], names(dgeResults_annotated[[i]])[1], "ENSEMBL")
    ens2symbol <-
      clusterProfiler::bitr(dgeResults_annotated[[i]]$ENSEMBL,
        fromType = "ENSEMBL",
        toType = "SYMBOL",
        OrgDb = "org.Mm.eg.db"
      )
    dgeResults_annotated[[i]] <- dplyr::full_join(dgeResults_annotated[[i]], ens2symbol)
  }
  return(dgeResults_annotated)
}
