#' count_matrix_loader
#'
#' @param file_type your raw data file type (presently optimized for txt)
#'
#' @return a count matrix containing raw count values for the experiment with
#' group annotations.

count_matrix_assembly <- function(file_type) {
  count_file <- fs::dir_ls(here::here("data-raw/"),
    regexp = file_type,
    recurse = TRUE
  )
  count_matrix <- openxlsx::read.xlsx(xlsxFile = count_file, sheet = 1)
  count_matrix$Geneid <- stringr::str_remove_all(count_matrix$Geneid, "\\..*")
  rownames(count_matrix) <- count_matrix$Geneid
  count_matrix <- count_matrix |>
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
    recurse = TRUE
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
    dplyr::filter(!Sample %in% outlier)
}

#' Remove identified outliers from counts and arrange in same order of metadata
#'
#' @param counts
#' @param metadata
#' @param outlier
#'
#' @return

removeOutlierCounts <- function(counts, metadata) {
  counts <- counts %>%
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

#' Prepare DGElist and filter by expression
#'
#' @param metadata
#' @param counts
#' @param design
#'
#' @return

DGEprep <- function(metadata, counts, design) {
  group <- as.matrix(metadata$Group)
  RNAseq <- edgeR::DGEList(
    counts = counts,
    group = group
  )
  keep <- edgeR::filterByExpr(RNAseq,
    design = design
  )
  RNAseq <- RNAseq[keep, , keep.lib.sizes = F]
  RNAseq <- edgeR::calcNormFactors(RNAseq)
  return(RNAseq)
}

#'  DGE analysis of counts using generated contrast matrix
#'
#' @param counts
#' @param metadata
#' @param design
#' @param ctrsts
#'
#' @return

DGEanalysis <- function(RNAseq, design, ctrsts) {
  RNAseq <- edgeR::estimateDisp(
    RNAseq,
    design
  )
  efit <- edgeR::glmQLFit(
    RNAseq,
    design
  )
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


#' MDS analysis of RNAseq data
#'
#' @param RNAseq
#'
#' @return

MDSanalysis <- function(RNAseq, metadata) {
  mdsData <- limma::plotMDS(RNAseq, plot = FALSE)
  varianceExplained <- mdsData$var.explained
  mdsData <-
    mdsData$eigen.vectors %>%
    as.data.table() %>%
    dplyr::mutate(ID = rownames(RNAseq$samples)) %>%
    dplyr::mutate(Group = metadata$Group) %>%
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
    ggplot2::ggtitle("MDS Plot - RNAseq analysis") +
    ggplot2::xlab(paste("Dim1 (", round(100 * varianceExplained[1], 2), " %)", sep = "")) +
    ggplot2::ylab(paste("Dim2 (", round(100 * varianceExplained[2], 2), " %)", sep = ""))
  return(pBase)
}

#' Generate upset plot with colors
#'
#' @param dgeResults_annotated
#'
#' @return An upset plot

UpsetplotGeneration <- function(dgeResults_annotated) {
  sig_genes_names <- names(dgeResults_annotated)
  sig_genes <- vector(mode = "list", length = 3)
  names(sig_genes) <- names(dgeResults_annotated)
  for (i in 1:3) {
    sig_genes[[i]] <- dgeResults_annotated[[i]]
    sig_genes[[i]] <- sig_genes[[i]] %>%
      dplyr::filter(FDR < 0.05) |>
      dplyr::mutate(
        Direction =
          dplyr::case_when(
            logFC > 0 ~ "Up",
            logFC < 0 ~ "Down"
          )
      )
    sig_genes[[i]] <- sig_genes[[i]] |>
      dplyr::select(SYMBOL, Direction)
  }


  names(sig_genes) <- c("Liver vs CS", "Liver vs PH", "CS vs PH")


  # make same upsetplot with ComplexUpset
  order_upset <- c("Liver vs CS", "Liver vs PH", "CS vs PH")
  upset_data <- data.frame(
    "Genes" = sig_genes[[1]]$SYMBOL,
    "Group" = names(sig_genes[1]),
    "Direction" = sig_genes[[1]]$Direction
  )

  for (i in 2:length(sig_genes)) {
    upset_data <- dplyr::add_row(upset_data,
      "Genes" = sig_genes[[i]]$SYMBOL,
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
     ggplot2::ggtitle("No. of differentially expressed genes between groups") +
    ggplot2::theme(plot.title = ggplot2::element_text(
      size = 24,
      hjust = 0.5,
      vjust = 0.95
    ))
  upsetRNA <- ggplotify::as.ggplot(upsetRNA)
  return(upsetRNA)
}

#' XY plot for CPM by group
#'
#' @param dgeList
#'
#' @return

XYCPMplot <- function(dgeList) {
  y <- edgeR::cpmByGroup(dgeList, log = T)
  cpm1 <- ggplot2::ggplot(
    as.data.frame(y),
    ggplot2::aes(x = L, y = PH)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::ggtitle("logCPM L vs PH") +
    ggplot2::theme(plot.title = ggplot2::element_text(
      hjust = 0.5,
      size = 18
    )) +
    ggplot2::xlab("logCPM L") +
    ggplot2::ylab("logCPM PH")

  cpm2 <- ggplot2::ggplot(
    as.data.frame(y),
    ggplot2::aes(x = CS, y = PH)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::ggtitle("logCPM CS vs PH") +
    ggplot2::theme(plot.title = ggplot2::element_text(
      hjust = 0.5,
      size = 18
    )) +
    ggplot2::xlab("logCPM CS") +
    ggplot2::ylab("logCPM PH")

  cpm_figs <- cpm1 + cpm2
  return(cpm_figs)
}

#' Helper function for GO analysis to subset data in significant genes,
#' split by increased or decreased expression
#'
#' @param dgeResults_annotated
#' @param Groupnames
#'
#' @return
#' @export
#'
#' @examples
PrepareComparison <- function(dgeResults_annotated, Groupnames) {
  dataList <- vector(mode = "list", length = 2)
  if ("FDR" %in% colnames(dgeResults_annotated)) {
    dataList[[1]] <- dgeResults_annotated |>
      dplyr::filter(logFC > 0) |>
      dplyr::filter(FDR < 0.05) |>
      dplyr::select(SYMBOL)
    dataList[[2]] <- dgeResults_annotated |>
      dplyr::filter(logFC < 0) |>
      dplyr::filter(FDR < 0.05) |>
      dplyr::select(SYMBOL)
  } else if ("adj.P.Val" %in% colnames(dgeResults_annotated)) {
    dataList[[1]] <- dgeResults_annotated |>
      dplyr::filter(logFC > 0) |>
      dplyr::filter(adj.P.Val < 0.05) |>
      dplyr::mutate(SYMBOL = Genes) |>
      dplyr::select(SYMBOL)
    dataList[[2]] <- dgeResults_annotated |>
      dplyr::filter(logFC < 0) |>
      dplyr::filter(adj.P.Val < 0.05) |>
      dplyr::mutate(SYMBOL = Genes) |>
      dplyr::select(SYMBOL)
  } else if ("padj" %in% colnames(dgeResults_annotated)) {
      dataList[[1]] <- dgeResults_annotated |>
          dplyr::filter(log2FoldChange > 0) |>
          dplyr::filter(padj < 0.05) |>
          dplyr::mutate(SYMBOL = gene) |>
          dplyr::select(SYMBOL)
      dataList[[2]] <- dgeResults_annotated |>
          dplyr::filter(log2FoldChange < 0) |>
          dplyr::filter(padj < 0.05) |>
          dplyr::mutate(SYMBOL = gene) |>
          dplyr::select(SYMBOL)
  }

  names(dataList) <- Groupnames
  dataList <- dataList |>
    purrr::map(~ clusterProfiler::bitr(
      .$SYMBOL,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = "org.Mm.eg.db",
      drop = T
    )[2]) |>
    purrr::map(~ .$ENTREZID)
  return(dataList)
}


#' GO analysis for cellular component, split by up- and down-regulated genes
#'
#' @param dgeResults_annotated
#'
#' @return

GOCCSplit <- function(dgeResults_annotated) {
  L_vs_CS <- PrepareComparison(
    dgeResults_annotated[[1]],
    c("Upregulated in L vs CS", "Upregulated in CS vs L")
  )
  L_vs_PH <- PrepareComparison(
    dgeResults_annotated[[2]],
    c("Upregulated in L vs PH", "Upregulated in PH vs L")
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
    ont = "CC",
    readable = T
  )
  return(clusterCompare_All)
}

#' Dotplot for CC GO analysis
#'
#' @param clusterCompare_All
#'
#' @return

DotplotCC <- function(clusterCompare_All, title) {
  CompareClusterFigure_All <- clusterProfiler::dotplot(clusterCompare_All,
                                                       by = "count",
                                                       showCategory = 5) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 30,
      vjust = 1,
      hjust = 1
    ),
    plot.title = ggplot2::element_text(size = 24,
                                       hjust = 0.5)) +
    ggplot2::xlab("") +
    ggplot2::ggtitle(title)
  return(CompareClusterFigure_All)
}

#' Upsetplot for GO terms
#'
#' @param clusterCompare_All
#'
#' @return

UpsetGO <- function(clusterCompare_All) {
  sig_GO <- clusterCompare_All@compareClusterResult |>
    dplyr::select(ID, Cluster)

  sig_GO <- sig_GO |>
    dplyr::mutate("TRUE" = TRUE)

  upset_GO_wide <- tidyr::pivot_wider(sig_GO, names_from = Cluster, values_from = "TRUE", values_fill = FALSE)

  order_upset_GO <- names(clusterCompare_All@geneClusters)

  upsetGO <- ComplexUpset::upset(upset_GO_wide, order_upset_GO, name = "", sort_sets = "descending", themes = ComplexUpset::upset_modify_themes(list(
    "intersections_matrix" = ggplot2::theme(text = ggplot2::element_text(size = 16))
  )))
  # ggplotify to use the object in patchwork
  upsetGO <- ggplotify::as.ggplot(upsetGO)
  upsetGO <- upsetGO +
    ggplot2::ggtitle("No. of GO terms from the \'cellular component\' ontology") +
    ggplot2::theme(plot.title = ggplot2::element_text(
      size = 24,
      hjust = 0.5,
      vjust = 0.95
    ))


  upsetGO
}
#' Function to create XY plots for CPM values for proteomics
#'
#' @param normalized_proteomics_res
#'
#' @return

CPMPlotsProteomics <- function(normalized_proteomics_res) {
  y <- as.data.frame(normalized_proteomics_res) |>
    dplyr::rowwise() |>
    dplyr::summarise(
      L = mean(c_across(Liver1:Liver8)),
      CS = mean(c_across(CS1:CS8)),
      PH = mean(c_across(PH1:PH8))
    ) |>
    dplyr::select(L, CS, PH)
  cpm1 <- ggplot2::ggplot(
    y,
    ggplot2::aes(x = L, y = PH)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::ggtitle("Abundance L vs PH") +
    ggplot2::theme(plot.title = ggplot2::element_text(
      hjust = 0.5,
      size = 18
    )) +
    ggplot2::xlab("Norm. Abun. L") +
    ggplot2::ylab("Norm. Abun. PH")

  cpm2 <- ggplot2::ggplot(
    y,
    ggplot2::aes(x = CS, y = PH)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::ggtitle("Abundance CS vs PH") +
    ggplot2::theme(plot.title = ggplot2::element_text(
      hjust = 0.5,
      size = 18
    )) +
    ggplot2::xlab("Norm. Abun. CS") +
    ggplot2::ylab("Norm. Abun. PH")

  cpm_figs <- cpm1 + cpm2
  return(cpm_figs)
}


#' Title
#'
#' @param GOobject
#' @param targetrow
#' @param DGElist
#' @param setup
#'
#' @return

RNAHeatmap <- function(GOobject,targetrow, DGElist, setup){
    #prepare setup file
    setup <-setup |>
            dplyr::arrange(Group = base::factor(Group, c("L", "CS", "PH")))
        #prepare CPM matrix from DEG list
        counts <- edgeR::cpm(DGElist,log = T)
        counts <- counts |>
            as.data.frame() |>
            dplyr::select(setup$Sample) |>
            dplyr::arrange(dplyr::desc(rowSums(counts)))
        counts <- counts |>
            dplyr::mutate(ENSEMBL = rownames(counts))

        #generate SYMBOL IDs for count matrix
        keyID<- clusterProfiler::bitr(counts$ENSEMBL,
                                    fromType = "ENSEMBL",
                                    toType = "SYMBOL",
                                    OrgDb = "org.Mm.eg.db",
                                    drop = T)
        keyID <- keyID |>
            dplyr::distinct(ENSEMBL, .keep_all = T)

        counts <- dplyr::left_join(counts,keyID, by = c("ENSEMBL"="ENSEMBL"))

        #Generate Gene list
        if(class(GOobject)=="compareClusterResult"){
            gene_list <- GOobject@compareClusterResult$geneID[targetrow]
            heatmapname <- GOobject@compareClusterResult$Description[targetrow]
        }
        else if(class(GOobject)=="enrichResult"){
            gene_list <- GOobject@result$geneID[targetrow]
            heatmapname <- GOobject@result$Description[targetrow]
        }
        else{
            print("No GO object detected")
        }

        gene_list <- unlist(stringr::str_split(gene_list, "/"))

        #Select Candidates in Gene List
        trimmed_cpm <- counts|>
            dplyr::select(-ENSEMBL) |>
            dplyr::filter(SYMBOL %in% gene_list) |>
            dplyr::distinct(SYMBOL, .keep_all = T) |>
            dplyr::arrange(SYMBOL)
        trimmed_cpm <- trimmed_cpm |>
            dplyr::filter(!is.na(SYMBOL))
        rownames(trimmed_cpm)<-trimmed_cpm$SYMBOL
        trimmed_cpm <- trimmed_cpm |>
            dplyr::select(-SYMBOL)

        #create annotation key for heatmap
        key <- as.data.frame(setup)
        key <- key |>
            dplyr::select(Group)
        rownames(key) <- setup$Sample
        key$Group <- factor(key$Group, c("L", "CS", "PH"))

        #create heatmap
        Heatmap_title <- grid::textGrob(heatmapname,
                                        gp = gpar(fontsize = 24,
                                                  fontface = "bold"))
        Heatmap <- pheatmap::pheatmap(trimmed_cpm,
                                      treeheight_col = 0,
                                      treeheight_row = 0,
                                      scale = "row",
                                      legend = T,
                                      na_col = "white",
                                      Colv = NA,
                                      na.rm = T,
                                      cluster_cols = F,
                                      fontsize_row = 5,
                                      fontsize_col = 8,
                                      cellwidth = 20,
                                      cellheight = 1.5,
                                      annotation_col = key,
                                      show_colnames = F,
                                      show_rownames = F,
                                      cluster_rows = T
        )
        Heatmap <- gridExtra::grid.arrange(grobs = list(Heatmap_title,
                                                        Heatmap[[4]]),
                                           heights = c(0.1, 1))
        Heatmap <- ggplotify::as.ggplot(Heatmap, scale = 1)
        Heatmap <- ggplotify::as.ggplot(Heatmap, scale = 1)
        return(Heatmap)



}
