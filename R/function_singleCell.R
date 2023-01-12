
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

#' Differential gene expression analysis on pseudocounts
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

#' Protein RNA correlation
#'
#' @param RNAdata
#' @param proteindata
#'
#' @return

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

#' Correlate proteins and genes with differential abundance and expression
#'
#' @param ProteinRNAComparison
#'
#' @return

ProteinRNACorFigure <- function(ProteinRNAComparison){
    ComparisonPlot <- ggplot2::ggplot(ProteinRNAComparison,
                                      ggplot2::aes(x = log2FoldChange,
                                                   y = logFC))+
        ggplot2::geom_point(size = 2)+
        ggplot2::theme_bw()+
        ggplot2::xlab("Pseudobulk RNA exp. Log2FC")+
        ggplot2::ylab("Protein Ab. Log2FC")+
        ggplot2::ggtitle("Liver vs PH",
                         "Protein abundance vs. RNA expression")+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                          size = 24),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5,
                                                             size = 22),
                       axis.title = ggplot2::element_text(size = 20),
                       axis.text = ggplot2::element_text(size = 20))
    return(ComparisonPlot)
}

#' Title
#'
#' @param ProtRNAComparison
#' @param limma_data the results file from the proteomics analysis to construct the background
#'
#' @return

GOCCRNAProt <- function(ProtRNAComparison, limma_data){
    bg <- clusterProfiler::bitr(limma_data$L_vs_CS$Genes,
                                fromType = "SYMBOL",
                                toType = "ENTREZID",
                                OrgDb = "org.Mm.eg.db"
    )

    sub_group_comparison <- list(Up_Prot_Up_RNA = dplyr::filter(ProtRNAComparison,
                                                                log2FoldChange > 0 &
                                                                    logFC > 0),
                                 Down_Prot_Down_RNA = dplyr::filter(ProtRNAComparison,
                                                                    log2FoldChange < 0 &
                                                                        logFC < 0),
                                 Down_Prot_Up_RNA = dplyr::filter(ProtRNAComparison,
                                                                  log2FoldChange > 0 &
                                                                      logFC < 0),
                                 Up_Prot_Down_RNA = dplyr::filter(ProtRNAComparison,
                                                                  log2FoldChange < 0 &
                                                                      logFC > 0))

    entrez_list <- lapply(sub_group_comparison,
                          function(x) clusterProfiler::bitr(x$gene,
                                                            fromType = "SYMBOL",
                                                            toType = "ENTREZID",
                                                            OrgDb = "org.Mm.eg.db") |>
                              dplyr::pull(ENTREZID))
    names(entrez_list)<- c("PH down Prot + RNA",
                           "PH up Prot + RNA",
                           "PH up Prot down RNA",
                           "PH down Prot up RNA")
    GO_results <- clusterProfiler::compareCluster(entrez_list,
                                                  fun = clusterProfiler::enrichGO,
                                                  keyType = "ENTREZID",
                                                  ont = "CC",
                                                  universe = bg$ENTREZID,
                                                  OrgDb = "org.Mm.eg.db",
                                                  readable = T)
    return(GO_results)
}

#' Plot MDS plot for dim1 and dim2 from DESeqDataset
#'
#' @param DESEQ2Object
#'
#' @return

MDSPseudoPlot <- function(DESEQCounts){
    mdsData <- limma::plotMDS(DESEQCounts@assays@data@listData[[1]],
                              plot = FALSE)
    varianceExplained <- mdsData$var.explained
    mdsData <-
        mdsData$eigen.vectors |>
        as.data.frame() |>
        dplyr::mutate(ID = DESEQCounts$hash.mcl.ID)  |>
        dplyr::mutate(Group = DESEQCounts$group)  |>
        dplyr::select(ID, Group, V1, V2, V3) |>
        dplyr::rename(dim1 = V1,
                      dim2 =  V2,
                      dim3 = V3)

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

#' Generate upset plot with colors from DESEQ objects
#'
#' @param dgeResults_annotated
#'
#' @return An upset plot

UpsetplotGenerationPseudo <- function(dgeResults_annotated) {
    sig_genes_names <- names(dgeResults_annotated)
    sig_genes <- vector(mode = "list", length = 3)
    names(sig_genes) <- names(dgeResults_annotated)
    for (i in 1:3) {
        sig_genes[[i]] <- dgeResults_annotated[[i]]
        sig_genes[[i]] <- sig_genes[[i]] %>%
            dplyr::filter(padj < 0.05) |>
            dplyr::mutate(
                Direction =
                    dplyr::case_when(
                        log2FoldChange > 0 ~ "Up",
                        log2FoldChange < 0 ~ "Down"
                    )
            )
        sig_genes[[i]] <- sig_genes[[i]] |>
            dplyr::select(gene, Direction)
    }


    names(sig_genes) <- c("Liver vs PH", "Liver vs CS", "CS vs PH")


    # make same upsetplot with ComplexUpset
    order_upset <- c("Liver vs CS", "Liver vs PH", "CS vs PH")
    upset_data <- data.frame(
        "Genes" = sig_genes[[1]]$gene,
        "Group" = names(sig_genes[1]),
        "Direction" = sig_genes[[1]]$Direction
    )

    for (i in 2:length(sig_genes)) {
        upset_data <- dplyr::add_row(upset_data,
                                     "Genes" = sig_genes[[i]]$gene,
                                     "Group" = names(sig_genes)[i],
                                     "Direction" = sig_genes[[i]]$Direction
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

    upsetRNA <- ComplexUpset::upset(upset_wide,
                                    order_upset,
                                    name = "",
                                    sort_sets = F,
                                    themes = ComplexUpset::upset_modify_themes(list(
                                        "intersections_matrix" = ggplot2::theme(text = ggplot2::element_text(size = 16))
                                    ))
    )
    # ggplotify to use the object in patchwork
    upsetRNA <- ggplotify::as.ggplot(upsetRNA)
    upsetRNA <- upsetRNA +
        # ggplot2::ggtitle("Upsetplot") +
        ggplot2::theme(plot.title = ggplot2::element_text(
            size = 18,
            hjust = 0.5,
            vjust = 0.95
        ))
    upsetRNA <- ggplotify::as.ggplot(upsetRNA)
    return(upsetRNA)
}

#' GO analysis for cellular component, split by up- and down-regulated genes
#'
#' @param dgeResults_annotated
#'
#' @return

GOCCSplitPseudo <- function(dgeResults_annotated) {
    L_vs_PH <- PrepareComparison(
        dgeResults_annotated[[1]],
        c("Upregulated in L vs PH", "Upregulated in PH vs L")
    )
    L_vs_CS <- PrepareComparison(
        dgeResults_annotated[[2]],
        c("Upregulated in L vs CS", "Upregulated in CS vs L")
    )
    PH_vs_CS <- PrepareComparison(
        dgeResults_annotated[[3]],
        c("Upregulated in CS vs PH", "Upregulated in PH vs CS")
    )
    All_comparisons <- c(L_vs_CS, L_vs_PH, PH_vs_CS)

    bg <- dgeResults_annotated[[1]]
    if ("Genes" %in% colnames(bg)) {
        bg <- bg |>
            dplyr::rename(SYMBOL = Genes)
    }
    else if ("gene" %in% colnames(bg)){
        bg <- bg |>
            dplyr::rename(SYMBOL = gene)
    }
    bg <- clusterProfiler::bitr(
        bg$SYMBOL,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        drop = T
    )
    clusterCompare_All <- clusterProfiler::compareCluster(
        gene = All_comparisons,
        fun = "enrichGO",
        universe = bg$ENTREZID,
        key = "ENTREZID",
        OrgDb = "org.Mm.eg.db",
        ont = "CC"
    )
    return(clusterCompare_All)
}
