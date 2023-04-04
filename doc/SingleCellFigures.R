#s-curve

#The UMAP plots are too large to save as target objects. Instead, they are
#generated and saved through this script
library(Brobdingnag)
#restart R after running Brobdingnag
library(tidyverse)
library(Seurat)
library(SummarizedExperiment)
source("~/R/tmp/PrimHepMultiOmics/R/functions.R")
source("~/R/tmp/PrimHepMultiOmics/R/functionsProteomics.R")
source("~/R/tmp/PrimHepMultiOmics/R/function_singleCell.R")
options(ggrepel.max.overlaps = Inf)

SingleCellData <- targets::tar_read(seuratObject)
Markers <- readRDS(here::here("data/markersUpdated.rds"))

#####Dotplot for cell markers####
Seurat::Idents(SingleCellData)<- "celltype"
marker_genes <- c("Alb",
                  "Cps1",
                  "Pck1",
                  "Ptprb",
                  "Stab2",
                  "Clec4f",
                  "Timd4",
                  "Nrxn1",
                  "Ank3",
                  "Ptprc",
                  "Dock2",
                  "Bicc1",
                  "Thsd4" )

DotplotFigure <- Seurat::DotPlot(SingleCellData,
                                 features = marker_genes)+
    ggplot2::ggtitle("Marker Genes for Cell Populations")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       vjust = 0.8,
                                                       size = 18),
                   axis.text.y = ggplot2::element_text(size = 18),
                   plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 24))


#####UMAP plot for group####
Seurat::Idents(SingleCellData) <- "Group"
DimPlotGroup <- Seurat::DimPlot(SingleCellData,
                        pt.size = 0.7,
                        label = T,
                        label.size = 6,
                        repel = T,
                        cols = wesanderson::wes_palette("Darjeeling1"),
                        label.box = T) +
    ggplot2::ggtitle("DimPlot By Group") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none"
    )

#####UMAP plot for cell type####
Seurat::Idents(SingleCellData)<- "celltype"
DimPlotCellType<- Seurat::DimPlot(SingleCellData,
                         pt.size = 0.7,
                         label = T,
                         label.size = 6,
                         repel = T,
                         cols = wesanderson::wes_palette(7,
                                                         name = "FantasticFox1",
                                                         type = "continuous"),
                         label.box = T) +
    ggplot2::ggtitle("DimPlot By Cell Type") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none"
    )
#####UMAP plot for individual####
Seurat::Idents(SingleCellData)<- "hash.mcl.ID"
SingleCellData <- Seurat::RenameIdents(SingleCellData,
                            "L1" = "L1",
                            "L2" = "L2",
                            "L3" = "L3",
                            "L4" = "L4",
                            "L5" = "L5",
                            "L6" = "L6",
                            "L7" = "L7",
                            "L8" = "L8",
                            "CS1" = "CS1",
                            "CD2" = "CS2",
                            "CS3" = "CS3",
                            "CS4" = "CS4",
                            "CS5" = "CS5",
                            "CS6" = "CS6",
                            "CS7" = "CS7",
                            "CS8" = "CS8",
                            "PH1" = "PH1",
                            "PH2" = "PH2",
                            "PH3" = "PH3",
                            "PH4" = "PH4",
                            "PH5" = "PH5",
                            "PH6" = "PH6",
                            "PH7" = "PH7",
                            "PH8" = "PH8"
)
DimPlotIndividual <- Seurat::DimPlot(SingleCellData,
                                 pt.size = 0.7,
                                 label = T,
                                 label.size = 6,
                                 repel = T,
                                 cols = wesanderson::wes_palette(24,
                                                                 name = "Zissou1",
                                                                 type = "continuous"),
                                 label.box = T) +
    ggplot2::ggtitle("DimPlot By individual") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none"
    )

#####Make Cyp2e1 and Hal plots####
Cyp2e1Plot <- Seurat::FeaturePlot(SingleCellData,
                                  features = "Cyp2e1",
                                  pt.size = 1)+
    ggplot2::ggtitle("Cyp2e1 - Central Marker")+
    ggplot2::theme(title = ggplot2::element_text(size = 22))

Cyp2f2Plot <- Seurat::FeaturePlot(SingleCellData,
                                  features = "Hal",
                                  pt.size = 1)+
    ggplot2::ggtitle("Hal - Portal Marker")+
    ggplot2::theme(title = ggplot2::element_text(size = 22))

#####Figure 4####

options(ggrepel.max.overlaps = Inf)

design_layout4 <- "
1133
2244
5566
####
"
patchworktest4 <- DimPlotGroup+
    DimPlotIndividual+
    DimPlotCellType+
    DotplotFigure+
    Cyp2e1Plot+
    Cyp2f2Plot+
    patchwork::plot_layout(design = design_layout4)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

     # grDevices::pdf(here::here("data/Figure4.pdf"), height = 20, width = 20)
     #  patchworktest4
     #  dev.off()

#####Hepatocyte Subset####
#Subset data to only include hepatocytes and cultured hepatocytes

# SeuratHepatocytes <- subset(SingleCellData,
#                             celltype == "Hepatocytes"|
#                                 celltype == "Cultured Cells")
# SeuratHepatocytes <- subset(SeuratHepatocytes, seurat_clusters != 12)

# SeuratHepatocytes <- Seurat::RunPCA(SeuratHepatocytes,
#                            npcs=35,
#                            seed.use=42,
#                            verbose=F)
#runUMAP
# SeuratHepatocytes  <- Seurat::RunUMAP(SeuratHepatocytes,
#                              dims = 1:35,
#                              seed.use = 43)
# SeuratHepatocytes  <- Seurat::FindNeighbors(SeuratHepatocytes,
#                                    reduction = "pca",
#                                    dims = 1:35)
#
# SeuratHepatocytes  <- Seurat::FindClusters(SeuratHepatocytes ,
#                                   random.seed = 42,
#                                   resolution = 1.2,
#                                   verbose = FALSE)

 # saveRDS(SeuratHepatocytes,
 #         here::here("data/HepatocytesSubset.rds"))

SeuratHepatocytes <- readRDS(here::here("data/HepatocytesSubset.rds"))

# hepatocytes_features <- Seurat::FindAllMarkers(SeuratHepatocytes,
#                                                only.pos = T,
#                                                min.pct = 0.25,
#                                                logfc.threshold = 0.25)

#saveRDS(hepatocytes_features, here::here("data/hepatocyte_markers.rds"))

#####Group plot for hepatocytes####

Idents(SeuratHepatocytes) <- "hash.mcl.ID"
SeuratHepatocytes <- RenameIdents(SeuratHepatocytes,
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
SeuratHepatocytes$Group <- Idents(SeuratHepatocytes)

DimplotGroupHepatocytes <- Seurat::DimPlot(SeuratHepatocytes,
                                           pt.size = 0.7,
                                           label = T,
                                           label.size = 6,
                                           repel = T,
                                           cols = wesanderson::wes_palette("Darjeeling1"),
                                           label.box = T) +
    ggplot2::ggtitle("DimPlot By Group \n for hepatocytes") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none"
    )

#####Dimplot for hepatocytes by individual####
Seurat::Idents(SeuratHepatocytes)<- "hash.mcl.ID"
SeuratHepatocytes <- Seurat::RenameIdents(SeuratHepatocytes,
                                       "L1" = "L1",
                                       "L2" = "L2",
                                       "L3" = "L3",
                                       "L4" = "L4",
                                       "L5" = "L5",
                                       "L6" = "L6",
                                       "L7" = "L7",
                                       "L8" = "L8",
                                       "CS1" = "CS1",
                                       "CD2" = "CS2",
                                       "CS3" = "CS3",
                                       "CS4" = "CS4",
                                       "CS5" = "CS5",
                                       "CS6" = "CS6",
                                       "CS7" = "CS7",
                                       "CS8" = "CS8",
                                       "PH1" = "PH1",
                                       "PH2" = "PH2",
                                       "PH3" = "PH3",
                                       "PH4" = "PH4",
                                       "PH5" = "PH5",
                                       "PH6" = "PH6",
                                       "PH7" = "PH7",
                                       "PH8" = "PH8"
)
DimPlotIndividualHepatocytes <- Seurat::DimPlot(SeuratHepatocytes,
                                     pt.size = 0.7,
                                     label = T,
                                     label.size = 6,
                                     repel = T,
                                     cols = wesanderson::wes_palette(24,
                                                                     name = "Zissou1",
                                                                     type = "continuous"),
                                     label.box = T) +
    ggplot2::ggtitle("DimPlot By individual \ for hepatocytes") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none"
    )

#####Zonation markers hepatocytes####
Cyp2e1Hepa <- Seurat::FeaturePlot(SeuratHepatocytes,
                                  features = "Cyp2e1",
                                  pt.size = 1)+
    ggplot2::ggtitle("Cyp2e1 - Central Marker")+
    ggplot2::theme(title = ggplot2::element_text(size = 22))

HalHepa <- Seurat::FeaturePlot(SeuratHepatocytes,
                                  features = "Hal",
                                  pt.size = 1)+
    ggplot2::ggtitle("Hal- Portal Marker")+
    ggplot2::theme(title = ggplot2::element_text(size = 22))

#####Figure 5####


design_layout5 <- "
1133
2244
"
patchworktest5 <-
    DimPlotIndividualHepatocytes+
    DimplotGroupHepatocytes+
    Cyp2e1Hepa+
    HalHepa+
    patchwork::plot_layout(design = design_layout5)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

   # grDevices::pdf(here::here("data/Figure5.pdf"), height = 20, width = 20)
   #  patchworktest5
   #  dev.off()

#####PseudoDEG on hepatocytes####
PseudoHepatocytes <- Pseudobulk(SeuratHepatocytes)
DEGHepatocytes <-DGEPseudo(PseudoHepatocytes)

UpsetHepatocytes <- UpsetplotGenerationPseudo(DEGHepatocytes,
                                              "Pseudobulk - Hepatocyte Subset")

GOCCPHvsLUp <- clusterProfiler::enrichGO(gene = subset(DEGHepatocytes$L_vs_PH,
                                                       log2FoldChange<0 & padj<0.05) |> dplyr::pull(gene),
                                         OrgDb = "org.Mm.eg.db",
                                         keyType = "SYMBOL",
                                         universe = DEGHepatocytes$L_vs_PH$gene,
                                         ont = "CC")

GOCCPHvsLUpTreeData <- enrichplot::pairwise_termsim(GOCCPHvsLUp)
GOCCPHvsLUpTreePlot <- enrichplot::treeplot(GOCCPHvsLUpTreeData, nWords = 0) +
    ggplot2::ggtitle("Clustering of genes with \n increased expression in PH vs L")+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 24,
                                                      hjust = 0.5))

GOCCLvsPHUp <- clusterProfiler::enrichGO(gene = subset(DEGHepatocytes$L_vs_PH,
                                                       log2FoldChange>0 & padj<0.05) |> dplyr::pull(gene),
                                         OrgDb = "org.Mm.eg.db",
                                         keyType = "SYMBOL",
                                         universe = DEGHepatocytes$L_vs_PH$gene,
                                         ont = "CC")

GOCCLvsPHUpTreeData <- enrichplot::pairwise_termsim(GOCCLvsPHUp)
GOCCLvsPHUpTreePlot <- enrichplot::treeplot(GOCCLvsPHUpTreeData,nWords = 0) +
    ggplot2::ggtitle("Clustering of genes with \n increased expression in L vs PH")+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 24,
                                                      hjust = 0.5))

GOCCLvsPHUpTreePlot + GOCCPHvsLUpTreePlot

#####Figure 6####


design_layout6 <- "
1111
2233
2233
2233
"
patchworktest6 <-
    UpsetHepatocytes+
    GOCCLvsPHUpTreePlot+
    GOCCPHvsLUpTreePlot+
    patchwork::plot_layout(design = design_layout6)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

  # grDevices::pdf(here::here("data/Figure6.pdf"), height = 20, width = 20)
  # patchworktest6
  # dev.off()

#####Analysis and clustering of PH group####
#  SeuratPH <- subset(SingleCellData,
#                              Group == "PH")
# SeuratPH  <- Seurat::RunPCA(SeuratPH ,
#                             npcs=35,
#                             seed.use=42,
#                             verbose=F)
#
# SeuratPH  <- Seurat::RunUMAP(SeuratPH ,
#                               dims = 1:35,
#                               seed.use = 43)
# SeuratPH  <- Seurat::FindNeighbors(SeuratPH ,
#                                     reduction = "pca",
#                                     dims = 1:35)
#
# SeuratPH   <- Seurat::FindClusters(SeuratPH  ,
#                                    random.seed = 42,
#                                    resolution = 1.2,
#                                    verbose = FALSE)

 #saveRDS(SeuratPH, here::here("data/PHSubset.rds"))
  SeuratPH <- readRDS(here::here("data/PHSubset.rds"))
 #
 #  PH_features <- Seurat::FindAllMarkers(SeuratPH ,
 #                                                 only.pos = T,
 #                                                 min.pct = 0.25,
 #                                                 logfc.threshold = 0.25)
 #
 # saveRDS(PH_features, here::here("data/PH_markers.rds"))
PH_markers <- readRDS(here::here("data/PH_markers.rds"))

DimplotClusterPH <- Seurat::DimPlot(SeuratPH,
                                           pt.size = 0.7,
                                           label = T,
                                           label.size = 6,
                                           repel = T,
                                           cols = wesanderson::wes_palette(22,
                                                                           name = "Darjeeling1",
                                                                           type = "continuous"),
                                           label.box = T) +
    ggplot2::ggtitle("DimPlot By Cluster \n for PH samples") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none"
    )

#assign cell identities based on marker gene expression
Seurat::Idents(SeuratPH)<- "seurat_clusters"
SeuratPH <- RenameIdents(SeuratPH,
                            "0" = "Dedifferentiated Hepatocytes",
                            "1" = "Dedifferentiated Hepatocytes",
                            "2" = "Dedifferentiated Hepatocytes",
                            "3" = "Dedifferentiated Hepatocytes",
                            "4" = "Dedifferentiated Hepatocytes",
                            "5" = "Dedifferentiated Hepatocytes",
                            "6" = "Dedifferentiated Endothelial Cells",
                            "7" = "Dedifferentiated Hepatocytes",
                            "8" = "Dedifferentiated Hepatocytes",
                            "9" = "Dedifferentiated Hepatocytes",
                            "10" = "Dedifferentiated Hepatocytes",
                            "11" = "Dedifferentiated Hepatocytes",
                            "12" = "Dedifferentiated Hepatocytes",
                            "13" = "Dedifferentiated Hepatocytes",
                            "14" = "Dedifferentiated Hepatocytes",
                            "15" = "Dedifferentiated Hepatocytes",
                            "16" = "Dedifferentiated Hepatocytes",
                            "17" = "Stellate Cells",
                            "18" = "Dedifferentiated Hepatocytes",
                            "19" = "Dedifferentiated Hepatocytes",
                            "20" = "Endothelial Cells",
                            "21" = "Hepatocytes"
)

SeuratPH$celltype <- Idents(SeuratPH)
marker_genes <- c("Alb",
                  "Cps1",
                  "Pck1",
                  "Ptprb",
                  "Stab2",
                  "Clec4f",
                  "Timd4",
                  "Nrxn1",
                  "Ank3",
                  "Arhgap31",
                  "Adgrf5",
                  "Serpine1",
                  "Gstp2")
DotplotFigurePH <- Seurat::DotPlot(SeuratPH,
                                 features = marker_genes)+
    ggplot2::ggtitle("Marker Genes - PH subset")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       vjust = 0.8,
                                                       size = 18),
                   axis.text.y = ggplot2::element_text(size = 18),
                   axis.title = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 24,
                                                      hjust = 0.5))

Idents(SeuratPH) <- "hash.mcl.ID"
DimplotIndividualPH <- Seurat::DimPlot(SeuratPH,
                pt.size = 0.7,
                label = T,
                label.size = 6,
                repel = T,
                cols = wesanderson::wes_palette(8,
                                                name = "Zissou1",
                                                type = "continuous"),
                label.box = T) +
    ggplot2::ggtitle("DimPlot By individual \ for PH") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none"
    )

Idents(SeuratPH) <- "celltype"
DimplotCelltypePH <- Seurat::DimPlot(SeuratPH,
                                  pt.size = 0.7,
                                  label = T,
                                  label.size = 6,
                                  repel = T,
                                  cols = wesanderson::wes_palette("Darjeeling1"),
                                  label.box = T) +
    ggplot2::ggtitle("DimPlot By cell type \n for PH") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none"
    )

#####Subset of endothelial cells####
#  SeuratEndothelial <- subset(SingleCellData,
#                              celltype == "Endothelial Cells"|
#                                  seurat_clusters == 12)
# SeuratEndothelial <- Seurat::RunPCA(SeuratEndothelial,
#                             npcs=35,
#                             seed.use=42,
#                             verbose=F)
# #runUMAP
# SeuratEndothelial  <- Seurat::RunUMAP(SeuratEndothelial,
#                               dims = 1:35,
#                               seed.use = 43)
# SeuratEndothelial  <- Seurat::FindNeighbors(SeuratEndothelial,
#                                     reduction = "pca",
#                                     dims = 1:35)
#
# SeuratEndothelial  <- Seurat::FindClusters(SeuratEndothelial ,
#                                    random.seed = 42,
#                                    resolution = 1.2,
#                                    verbose = FALSE)

 # saveRDS(SeuratEndothelial,
 #         here::here("data/EndothelialSubset.rds"))

 SeuratEndothelial<- readRDS(here::here("data/EndothelialSubset.rds"))
 Idents(SeuratEndothelial)<- "Group"

 DimplotEndothelial <- Seurat::DimPlot(SeuratEndothelial,
                                        pt.size = 0.7,
                                        label = T,
                                        label.size = 6,
                                        repel = T,
                                        cols = wesanderson::wes_palette(3,
                                                                        name = "Zissou1",
                                                                        type = "continuous"),
                                        label.box = T) +
     ggplot2::ggtitle("DimPlot for Endothelial Cells") +
     ggplot2::theme(
         plot.title = ggplot2::element_text(hjust = 0.5,
                                            size = 24),
         legend.position = "none"
     )
DimplotEndothelial
Idents(SeuratEndothelial)<- "seurat_clusters"
# Endothelial_features <- Seurat::FindAllMarkers(SeuratEndothelial,
#                                       only.pos = T,
#                                       min.pct = 0.25,
#                                       logfc.threshold = 0.25)
#
# saveRDS(Endothelial_features, here::here("data/EndothelialFeatures.rds"))

Endothelial_features <- readRDS(here::here("data/EndothelialFeatures.rds"))

markers_endofigure <- c("Arhgap31",
                        "Adgrf5",
                        "Gata4",
                        "Itgb3",
                        "Pkm",
                        "Gclc",
                        "Lipc")
Idents(SeuratEndothelial)<- "Group"

DotplotFigureEndo <- Seurat::DotPlot(SeuratEndothelial,
                                 features = markers_endofigure)+
    ggplot2::ggtitle("Marker Genes - Endothelial cells")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       vjust = 0.8,
                                                       size = 18),
                   axis.text.y = ggplot2::element_text(size = 18),
                   plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 24),
                   axis.title = ggplot2::element_blank())


#Annotate dedifferentiated epithelial cell cluster
Seurat::Idents(SingleCellData)<-"seurat_clusters"
SingleCellData <- RenameIdents(SingleCellData,
                            "0" = "Dedifferentiated Hepatocyte",
                            "1" = "Dedifferentiated Hepatocyte",
                            "2" = "Hepatocytes",
                            "3" = "Endothelial Cells",
                            "4" = "Hepatocytes",
                            "5" = "Hepatocytes",
                            "6" = "Hepatocytes",
                            "7" = "Hepatocytes",
                            "8" = "Dedifferentiated Hepatocyte",
                            "9" = "Hepatocytes",
                            "10" = "Hepatocytes",
                            "11" = "Hepatocytes",
                            "12" = "Dedifferentiated Endothelial Cells",
                            "13" = "Hepatocytes",
                            "14" = "Dedifferentiated Hepatocyte",
                            "15" = "Hepatocytes",
                            "16" = "Kupfer Cells",
                            "17" = "Dedifferentiated Hepatocyte",
                            "18" = "Hepatocytes",
                            "19" = "Hepatocytes",
                            "20" = "Stellate Cells",
                            "21" = "Hepatocytes",
                            "22" = "Hepatocytes",
                            "23" = "Hepatocytes",
                            "24" = "Endothelial Cells",
                            "25" = "Leukocyte",
                            "26" = "Endothelial Cells",
                            "27" = "Stellate Cells",
                            "28" = "Billiary Epithelial Cells",
                            "29" = "Leukocyte",
                            "30" = "Billiary Epithelial Cells"
)
SingleCellData$updatedcelltype <- Seurat::Idents(SingleCellData)
cellcount <- SingleCellData@meta.data |>
    dplyr::group_by(Group) |>
    dplyr::select(Group, updatedcelltype) |>
    dplyr::count(updatedcelltype) |>
    dplyr::rename(Celltype = updatedcelltype)

cellcountfigure <- ggplot2::ggplot(cellcount, ggplot2::aes(x = n,
                                        y = Group,
                                        fill = Celltype)) +
    ggplot2::geom_bar(position = "stack",
                      stat = "identity") +
    ggplot2::ggtitle("Cell Distribution in Sample Types") +
    ggplot2::theme(
        title = ggplot2::element_text(size = 24,
                                      hjust = 0.5),
        axis.text.y = ggplot2::element_text(size = 16),
        axis.text.x = ggplot2::element_text(size = 14),
        legend.text = ggplot2::element_text(size = 12),
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "black")

    )+
    ggplot2::ylab("")+
    ggplot2::xlab("No. of cells")+
    ggplot2::scale_fill_manual(
        values = wesanderson::wes_palette("BottleRocket2",
                                          n = 8,
                                          type = "continuous")
    )
#make same plot but with ratio
cellcount <- cellcount |>
    dplyr::group_by(Group) |>
    dplyr::mutate(
        totalcells = sum(n),
        percentage = n/totalcells*100
    )
cellratio <- ggplot2::ggplot(cellcount, ggplot2::aes(x = percentage,
                                                     y = Group,
                                                     fill = Celltype)) +
    ggplot2::geom_bar(position = "stack",
                      stat = "identity") +
    ggplot2::ggtitle("% cells in Sample Types") +
    ggplot2::theme(
        title = ggplot2::element_text(size = 24,
                                      hjust = 0.5),
        axis.text.y = ggplot2::element_text(size = 16),
        axis.text.x = ggplot2::element_text(size = 14),
        legend.text = ggplot2::element_text(size = 12),
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "black")
    )+
    ggplot2::ylab("")+
    ggplot2::xlab("%cells")+
    ggplot2::scale_fill_manual(
        values = wesanderson::wes_palette("Darjeeling1",
                                          n = 8,
                                          type = "continuous"))


#####Figure 7####
design_layout7 <- "
1122
3344
5566
"
patchworktest7 <-
    DimplotCelltypePH+
    DotplotFigurePH+
    cellcountfigure+
    cellratio+
    DimplotEndothelial+
    DotplotFigureEndo+
    patchwork::plot_layout(design = design_layout7)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

     # grDevices::pdf(here::here("data/Figure7.pdf"), height = 20, width = 20)
     # patchworktest7
     # dev.off()

#####Protein and RNA correlation####
correlationData <- targets::tar_read(ProtRNACor)
clusterCelltypeKey <- SingleCellData@meta.data |>
    dplyr::select(celltype, seurat_clusters) |>
    dplyr::distinct()
annotatedMarkers <- dplyr::left_join(Markers,
                                     clusterCelltypeKey,
                                     by = c("cluster" = "seurat_clusters"))
annotatedCor <- dplyr::left_join(correlationData,
                                 annotatedMarkers,
                                 by = "gene") |>
    dplyr::select(gene, log2FoldChange, logFC, celltype) |>
    dplyr::distinct()

#select candidates with either high or low abundance in cultured hepatocytes
candidates <- annotatedCor|>
    dplyr::filter(logFC>2|logFC< -2) |>
    dplyr::arrange(logFC)

#make heatmap of candidates
proteomicsData <- targets::tar_read(NormalizedMatrix)
nameConverter <- clusterProfiler::bitr(rownames(proteomicsData),
                                       fromType = "ACCNUM",
                                       toType = "SYMBOL",
                                       OrgDb = "org.Mm.eg.db")
proteomicsDataAnno <- proteomicsData |>
    as.data.frame() |>
    dplyr::mutate(ACCID = rownames(proteomicsData))
proteomicsDataAnno <- dplyr::left_join(proteomicsDataAnno,
                                       nameConverter,
                                       by = c("ACCID"="ACCNUM")) |>
    dplyr::filter(!is.na(SYMBOL)) |>
    dplyr::distinct(SYMBOL, .keep_all = T)

proteomicsDataAnno <- proteomicsDataAnno|>
    magrittr::set_rownames(proteomicsDataAnno$SYMBOL) |>
    dplyr::select(-SYMBOL,
                  -ACCID)

#Select Candidates in Gene List
trimmed_cpm <- dplyr::filter(proteomicsDataAnno,
                             rownames(proteomicsDataAnno) %in% candidates$gene)

#create annotation key for heatmap
setup <- targets::tar_read(ProteomicsSetup)
key <- as.data.frame(setup)
key <- key |>
    dplyr::select(Tissue)
rownames(key) <- setup$SampleID
key$Tissue <- factor(key$Tissue, c("liver", "CS", "PH"))

#create heatmap
Heatmap_title <- grid::textGrob("Candidates from proteomics data",
                                gp = grid::gpar(fontsize = 24,
                                          fontface = "bold"))
HeatmapProteome <- pheatmap::pheatmap(trimmed_cpm,
                              treeheight_col = 0,
                              treeheight_row = 0,
                              scale = "row",
                              legend = T,
                              na_col = "white",
                              Colv = NA,
                              na.rm = T,
                              cluster_cols = F,
                              fontsize_row = 12,
                              fontsize_col = 12,
                              cellwidth = 16,
                              cellheight = 12,
                              annotation_col = key,
                              show_colnames = F,
                              show_rownames = T,
                              cluster_rows = T
)
HeatmapProteome <- gridExtra::grid.arrange(grobs = list(Heatmap_title,
                                                HeatmapProteome[[4]]),
                                   heights = c(0.1, 1))
HeatmapProteome <- ggplotify::as.ggplot(HeatmapProteome,scale = 1)

glut1Plot <- Seurat::FeaturePlot(SingleCellData,
                         features = "Slc2a1",
                         pt.size = 1)+
    ggplot2::ggtitle("Glut1 Expression\n Cultured Hepatocytes")+
    ggplot2::theme(title = ggplot2::element_text(size = 22))

glut2Plot <- Seurat::FeaturePlot(SingleCellData,
                                 features = "Slc2a2",
                                 pt.size = 1)+
    ggplot2::ggtitle("Glut2 Expression\n Freshly Isolated Hepatocytes")+
    ggplot2::theme(title = ggplot2::element_text(size = 22))

Idents(SingleCellData)<-"celltype"
DotplotCandidateFigure <- Seurat::DotPlot(SingleCellData,
                                 features = rownames(trimmed_cpm))+
    ggplot2::ggtitle("Candidate Genes - RNA expression")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       vjust = 0.8,
                                                       size = 18),
                   axis.text.y = ggplot2::element_text(size = 18),
                   plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 24)
                   )+
    ggplot2::coord_flip()


#####Figure 8####
design_layout8 <- "
1122
3344
"
patchworktest8 <-
    HeatmapProteome+
    DotplotCandidateFigure+
    glut1Plot+
    glut2Plot+
    patchwork::plot_layout(design = design_layout8)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

     grDevices::pdf(here::here("data/Figure8.pdf"), height = 20, width = 20)
     patchworktest8
     dev.off()

#####Figure 9####
#Figure 9 is a GO analysis of the terms that do not change in hepatocytes.
#tested on proteomics CC and RNA DEGs but they had few terms changing. Ended
#using MF for proteomics instead
DAPResults <- targets::tar_read(DAPResults)
nonsigprots <- DAPResults$L_vs_PH |>
    dplyr::filter(adj.P.Val > 0.05)

eg <- clusterProfiler::bitr(
    nonsigprots$Genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

bg <- clusterProfiler::bitr(
    DAPResults$L_vs_PH$Genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
)

unchangedProtMF <- clusterProfiler::enrichGO(
    gene = eg$ENTREZID,
    universe = bg$ENTREZID,
    key = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    ont = "MF",
    readable = T
)

UnchangedDotPlot <- clusterProfiler::dotplot(unchangedProtMF)+
    ggplot2::theme(
        plot.title = ggplot2::element_text(size = 24,
                                   hjust = 0.5)) +
    ggplot2::ggtitle("Proteins with unchanged abundance between L and PH")
NormalizedMatrix <- targets::tar_read(NormalizedMatrix)
ProteomicsSetup <- targets::tar_read(ProteomicsSetup)
CatActHm <- proteomicsHeatmap(unchangedProtMF,
                  targetrow = 1,
                  counts = NormalizedMatrix,
                  ProteomicsSetup,
                  2.5)

patchworktest9 <-
    (UnchangedDotPlot+patchwork::plot_spacer())/
    (CatActHm+patchwork::plot_spacer())+
    patchwork::plot_annotation(tag_levels = "A")

   # grDevices::pdf(here::here("data/Figure9.pdf"), height = 20, width = 20)
   # patchworktest9
   # dev.off()

