#' Import raw data from quant
#'
#' @param filename
#'
#' @return

ProteomicsDataLoader <- function(filename) {
  Proteomics_file <- paste(here::here("data-raw/"), filename, sep = "")

  Proteomics_dataset <- vroom::vroom(Proteomics_file,
    col_types = vroom::cols(
      Liver1 = vroom::col_double(),
      Liver2 = vroom::col_double(),
      Liver3 = vroom::col_double(),
      Liver4 = vroom::col_double(),
      Liver5 = vroom::col_double(),
      Liver6 = vroom::col_double(),
      Liver7 = vroom::col_double(),
      Liver8 = vroom::col_double(),
      CS1 = vroom::col_double(),
      CS2 = vroom::col_double(),
      CS3 = vroom::col_double(),
      CS4 = vroom::col_double(),
      CS5 = vroom::col_double(),
      CS6 = vroom::col_double(),
      CS7 = vroom::col_double(),
      CS8 = vroom::col_double(),
      PH1 = vroom::col_double(),
      PH2 = vroom::col_double(),
      PH3 = vroom::col_double(),
      PH4 = vroom::col_double(),
      PH5 = vroom::col_double(),
      PH6 = vroom::col_double(),
      PH7 = vroom::col_double(),
      PH8 = vroom::col_double(),
      Protein.Group = vroom::col_character(),
      Protein.Ids = vroom::col_character(),
      Protein.Names = vroom::col_character(),
      Genes = vroom::col_character(),
      First.Protein.Description = vroom::col_character()
    ),
    na = "NaN"
  )

  # remove first row
  Proteomics_dataset <- Proteomics_dataset[-1, ]
  # create data matrix. Find rowsum to keep duplicates with highest abundance
  Proteomics_processed <- Proteomics_dataset

  Proteomics_processed <- Proteomics_processed |>
    dplyr::rowwise() |>
    dplyr::mutate(rowsum = sum(c_across(Liver1:PH8), na.rm = T))

  Proteomics_processed <- Proteomics_processed |>
    dplyr::arrange(dplyr::desc(rowsum)) |>
    dplyr::distinct(Protein.Ids, .keep_all = T)
  return(Proteomics_processed)
}

IDkey <- function(count_matrix) {
  convKey <- count_matrix |>
    dplyr::select(Protein.Ids, Genes)
  return(convKey)
}

#' Generate a setup file for proteomics
#'
#' @return
#' @export
#'
#' @examples
GenerateSetup <- function(Proteomics_dataset) {
  setup <- data.frame(
    "SampleID" = colnames(Proteomics_dataset)[1:24],
    "Tissue" = NA
  )
  setup$Tissue[1:8] <- "liver"
  setup$Tissue[9:16] <- "CS"
  setup$Tissue[17:24] <- "PH"
  return(setup)
}

#' Normalize matrix and remove lowly abundant proteins
#'
#' @param Proteomics_matrix
#' @param setup
#'
#' @return

DataFiltering <- function(Proteomics_matrix,
                          setup) {
  Proteomics_processed <- Proteomics_matrix
  Proteomics_matrix <- as.matrix(Proteomics_matrix[, 1:24])
  rownames(Proteomics_matrix) <- Proteomics_processed$Protein.Ids
  normalized_proteomics <- limma::normalizeBetweenArrays(
    log(Proteomics_matrix),
    method = "quantile"
  )
  # Remove samples where more than 4 are missing pr group
  missingSamples_proteomics <- data.table::data.table(
    is.na(normalized_proteomics),
    keep.rownames = TRUE
  ) |>
    data.table::melt(
      measure.vars = colnames(normalized_proteomics),
      variable.name = "SampleID"
    )
  missingSamples_proteomics <- merge(setup,
    missingSamples_proteomics,
    by = "SampleID"
  )
  data.table::setnames(missingSamples_proteomics, "rn", "Accession")
  missingSamples_proteomics <- missingSamples_proteomics |>
    dplyr::group_by(Accession, Tissue) |>
    dplyr::mutate(nMissing = sum(value)) |>
    reshape2::dcast(Accession ~ Tissue,
      value.var = "nMissing",
      fun.aggregate = mean
    )

  cutoff <- 4
  tooManyMissing_proteomics <- missingSamples_proteomics %>%
    dplyr::filter(liver > cutoff |
      CS > cutoff |
      PH > cutoff)

  normalized_proteomics_res <- normalized_proteomics[!(rownames(normalized_proteomics) %in% tooManyMissing_proteomics$Accession), ] # nolint
  return(normalized_proteomics_res)
}

#' Run differentially abundant protein analysis using limma
#'
#' @param normalized_proteomics_res
#' @param setup
#'
#' @return

LimmaAnalysis <- function(normalized_proteomics_res, setup, key) {
  design <- stats::model.matrix(~ 0 + Tissue, setup)
  colnames(design) <- stringr::str_remove_all(colnames(design), "Tissue")
  fit <- limma::lmFit(normalized_proteomics_res, design = design, method = "robust")
  cont.matrix <- limma::makeContrasts(
    L_vs_CS = liver - CS,
    L_vs_PH = liver - PH,
    CS_vs_PH = CS - PH,
    levels = design
  )
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

  resultTables <- list(
    L_vs_CS = limma::topTable(fit2, coef = "L_vs_CS", number = Inf, p.value = 1) %>% data.table::data.table(keep.rownames = TRUE),
    L_vs_PH = limma::topTable(fit2, coef = "L_vs_PH", number = Inf, p.value = 1) %>% data.table::data.table(keep.rownames = TRUE),
    CS_vs_PH = limma::topTable(fit2, coef = "CS_vs_PH", number = Inf, p.value = 1) %>% data.table::data.table(keep.rownames = TRUE)
  )

  lapply(resultTables, data.table::setnames, "rn", "Protein.Ids")
  resultTables <- resultTables |>
    purrr::map(~ dplyr::left_join(., key, by = c("Protein.Ids" = "Protein.Ids")))
}
#' Make MDS plot for proteomics analysis
#'
#' @param normalized_proteomics_res
#' @param setup
#'
#' @return

ProteomicsMDS <- function(normalized_proteomics_res, setup) {
  mdsData <- limma::plotMDS(normalized_proteomics_res, plot = FALSE)
  varianceExplained <- mdsData$var.explained
  mdsData <-
    mdsData$eigen.vectors %>%
    as.data.table() %>%
    dplyr::mutate(ID = setup$SampleID) %>%
    dplyr::mutate(Group = setup$Tissue) %>%
    dplyr::select(ID, Group, V1, V2, V3)

  setnames(
    mdsData,
    c("V1", "V2", "V3", "ID", "Group"),
    c("dim1", "dim2", "dim3", "ID", "Group")
  )

  pBase <-
    ggplot2::ggplot(mdsData, ggplot2::aes(x = dim1, y = dim2, colour = Group)) +
    ggplot2::geom_point(size = 10) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 18),
      plot.title = ggplot2::element_text(size = 22, hjust = 0.5)
    ) +
    ggplot2::ggtitle("MDS Plot") +
    ggplot2::xlab(paste("Dim1 (", round(100 * varianceExplained[1], 2), " %)", sep = "")) +
    ggplot2::ylab(paste("Dim2 (", round(100 * varianceExplained[2], 2), " %)", sep = ""))
  return(pBase)
}


##### GO analysis####
# GO analysis function is present in the functions.r script

UpsetProteomics <- function(limma_data) {
  for (i in 1:length(limma_data)) {
    limma_data[[i]] <- limma_data[[i]] |>
      dplyr::filter(adj.P.Val < 0.05) |>
      dplyr::select(Genes)
    limma_data[[i]] <- limma_data[[i]]$Genes
  }
  names(limma_data) <- c("Liver vs CS",
                         "Liver vs PH",
                         "CS vs PH")
  order_upset <- c("Liver vs CS", "Liver vs PH", "CS vs PH")
  upset_data <- data.frame(
    "Genes" = limma_data[[1]],
    "Group" = names(limma_data[1])
  )

  for (i in 2:length(limma_data)) {
    upset_data <- dplyr::add_row(upset_data,
      "Genes" = limma_data[[i]],
      "Group" = names(limma_data)[i]
    )
  }
  upset_data <- dplyr::distinct(upset_data)
  upset_data <- upset_data |>
    dplyr::mutate("TRUE" = TRUE)
  upset_wide <- tidyr::pivot_wider(upset_data,
                                   names_from = Group,
                                   values_from = "TRUE",
                                   values_fill = FALSE) |>
    dplyr::filter(!is.na(Genes))

  rownames(upset_wide) <- upset_wide$Genes

  upset_wide <- dplyr::select(upset_wide, -Genes)

  upsetProt <- ComplexUpset::upset(upset_wide,
    order_upset,
    name = "",
    sort_sets = F,
    themes = ComplexUpset::upset_modify_themes(list(
      "intersections_matrix" = ggplot2::theme(text = ggplot2::element_text(size = 16))
    ))
  )
  # ggplotify to use the object in patchwork
  upsetProt <- ggplotify::as.ggplot(upsetProt)

  return(upsetProt)
}
