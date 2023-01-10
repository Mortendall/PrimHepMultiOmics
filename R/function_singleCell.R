
#' Import Single Cell data
#'
#' @param seuratobject prepared by SCOP
#'
#' @return

loadSingleCell <- function(seuratobject) {
  # read in data file
  data_file <- fs::dir_ls(here::here("data-raw/"),
    regexp = seuratobject,
    recurse = TRUE
  )
  seurat_test <- readRDS(data_file)

  # Assign samples to groups
  seurat_test <- readRDS(data_file)
  Idents(seurat_test) <- "hash.mcl.ID"
  seurat_test <- RenameIdents(seurat_test,
    "L1" = "L",
    "L2" = "L",
    "L3" = "L",
    "L4" = "L",
    "L5" = "L",
    "L6" = "L",
    "L7" = "L",
    "L8" = "L",
    "CS1" = "CS",
    "CD2" = "CS",
    "CS3" = "CS",
    "CS4" = "CS",
    "CS5" = "CS",
    "CS6" = "CS",
    "CS7" = "CS",
    "CS8" = "CS",
    "PH1" = "PH",
    "PH2" = "PH",
    "PH3" = "PH",
    "PH4" = "PH",
    "PH5" = "PH",
    "PH6" = "PH",
    "PH7" = "PH",
    "PH8" = "PH"
  )
  seurat_test$Group <- Idents(seurat_test)
  Idents(seurat_test) <- "seurat_clusters"

  # running the find markers function can be quite time consuming. Instead, I run
  # it once outside of the script and save the result.

  # Markers <- Seurat::FindAllMarkers(seurat_object,
  #                                   only.pos = T,
  #                                   min.pct = 0.25,
  #                                   logfc.threshold = 0.25)
  # saveRDS(Markers, here::here(data/markersUpdated.rds))

  # assign identities based on marker gene expression
  seurat_test <- RenameIdents(seurat_test,
    "0" = "Cultured Cells",
    "1" = "Cultured Cells",
    "2" = "Hepatocytes",
    "3" = "Endothelial Cells",
    "4" = "Hepatocytes",
    "5" = "Hepatocytes",
    "6" = "Hepatocytes",
    "7" = "Hepatocytes",
    "8" = "Cultured Cells",
    "9" = "Hepatocytes",
    "10" = "Hepatocytes",
    "11" = "Hepatocytes",
    "12" = "Cultured Cells",
    "13" = "Hepatocytes",
    "14" = "Cultured Cells",
    "15" = "Hepatocytes",
    "16" = "Kupfer Cells",
    "17" = "Cultured Cells",
    "18" = "Hepatocytes",
    "19" = "Hepatocytes",
    "20" = "Stellate Cells",
    "21" = "Hepatocytes",
    "22" = "Hepatocytes",
    "23" = "Hepatocytes",
    "24" = "Endothelial Cells",
    "25" = "Leukocyte",
    "26" = "Endothelial Cells",
    "27" = "Cultured Cells",
    "28" = "Billiary Epithelial Cells",
    "29" = "Leukocyte",
    "30" = "Billiary Epithelial Cells"
  )
  seurat_test$celltype <- Seurat::Idents(seurat_test)

  return(seurat_test)
}

#' Convert single-nucleus RNAseq to pseudobulk
#'
#' @param seuratObject
#'
#' @return a Deseqobject

Pseudobulk <- function(seuratObject) {
  Seurat::Idents(seuratObject) <- "Group"
  counts <- seuratObject@assays$RNA@counts
  metadata <- seuratObject@meta.data
  metadata$Group <- factor(metadata$Group)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = metadata
  )
  groups <- colData(sce)[, c("Group", "hash.mcl.ID")]
  # named vector of cluster id

  kids <- purrr::set_names(levels(sce$Group))

  # named vector of sample names
  sce$hash.mcl.ID <- factor(sce$hash.mcl.ID)
  sids <- purrr::set_names(levels(sce$hash.mcl.ID))

  # no of samples
  ns <- length(sids)

  # determine no. of cells pr sample
  n_cells <- as.numeric(table(sce$hash.mcl.ID))

  ## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
  m <- match(sids, sce$hash.mcl.ID)

  ## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
  ei <- data.frame(colData(sce)[m, ],
    n_cells,
    row.names = NULL
  ) |>
    select(-celltype, -seurat_clusters)

  # perform QC for sanity
  sce <- scuttle::addPerCellQC(sce)
  # identify outliers
  sce$is_outlier <- scuttle::isOutlier(metric = sce$total, nmads = 2, type = "both", log = T)
  sce <- sce[, !sce$is_outlier]
  # remove genes with fewer than 10 counts
  sce <- sce[rowSums(counts(sce) > 1) >= 10]

  # aggregate counts pr hash.id and cluster id

  groups <- colData(sce)[, c("Group", "hash.mcl.ID")]

  pb <- Matrix.utils::aggregate.Matrix(t(SingleCellExperiment::counts(sce)),
    groupings = groups,
    fun = "sum"
  )

  # Not every cluster is present in all samples; create a vector that represents how to split samples
  splitf <- sapply(
    stringr::str_split(rownames(pb),
      pattern = "_",
      n = 2
    ),
    `[`, 1
  )

  # Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
  pb <- magrittr::set_colnames(
    t(pb),
    stringr::str_extract(rownames(pb), "(?<=_)[:alnum:]+")
  )
  analysis_metadata <- data.frame(
    group = ei$Group,
    hash.mcl.ID = ei$hash.mcl.ID
  )


  rownames(analysis_metadata) <- analysis_metadata$hash.mcl.ID

  cluster_counts <- data.frame(pb)

  all(rownames(analysis_metadata) == colnames(cluster_counts))

  # create dds object
  dds <- DESeq2::DESeqDataSetFromMatrix(cluster_counts,
    colData = analysis_metadata,
    design = ~group
  )
  dds <- DESeq2::DESeq(dds)
  return(dds)
}

#' Create a count matrix from DESeq object
#'
#' @param DGEObject
#'
#' @return

DGECountMatrix <- function(DGEObject) {
  rld <- DESeq2::rlog(DGEObject,
    blind = T
  )
  return(rld)
}

#' Title
#'
#' @param DGEObject
#'
#' @return a differential gene expression analysis

DGEPseudo <- function(DGEObject) {
  contrast_list <- list(
    L_vs_PH = c(
      "group",
      levels(DGEObject$group)[1],
      levels(DGEObject$group)[3]
    ),
    L_vs_CS = c(
      "group",
      levels(DGEObject$group)[1],
      levels(DGEObject$group)[2]
    ),
    CS_vs_PH = c(
      "group",
      levels(DGEObject$group)[2],
      levels(DGEObject$group)[3]
    )
  )
  resList <- lapply(contrast_list, function(x) {
    DESeq2::results(DGEObject,
      contrast = x,
      alpha = 0.05
    )
  })
  res_tbl_list <- lapply(resList, function(x) {
    data.frame(x) |>
      tibble::rownames_to_column(var = "gene") |>
      dplyr::arrange(padj) |>
      tibble::as_tibble()
  })
  return(res_tbl_list)
}

ProteinRNACorrelation <- function(RNAdata, proteindata) {
  for (i in 1:length(proteindata)) {
    proteindata[[i]] <- proteindata[[i]] |>
      dplyr::filter(adj.P.Val < 0.05)
  }

  proteomics_L_vs_PH <- proteindata$L_vs_PH

  RNAseqPseudo_L_vs_PH <- RNAdata$L_vs_PH |>
    dplyr::filter(padj < 0.05)
  RNAseqPseudo_L_vs_PH_overlap <- RNAseqPseudo_L_vs_PH |>
    dplyr::filter(gene %in% proteomics_L_vs_PH$Genes)
  proteomics_L_vs_PH_overlap <- proteomics_L_vs_PH |>
    dplyr::filter(Genes %in% RNAseqPseudo_L_vs_PH_overlap$gene)

  joined_table <- dplyr::left_join(RNAseqPseudo_L_vs_PH_overlap,
    proteomics_L_vs_PH_overlap,
    by = c("gene" = "Genes")
  )
  return(joined_table)
}
