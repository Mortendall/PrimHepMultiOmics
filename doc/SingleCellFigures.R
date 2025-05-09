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
                         cols = wesanderson::wes_palette(8,
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
SingleCellData$hash.mcl.ID <- Idents(SingleCellData)
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
    ggplot2::ylab("")+
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none",
    )

#####Make Cyp2e1 and Hal plots####
Cyp2e1Plot <- Seurat::FeaturePlot(SingleCellData,
                                  features = "Cyp2e1",
                                  pt.size = 1)+
    ggplot2::ggtitle("Cyp2e1 - Central Marker")+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 22))

Cyp2f2Plot <- Seurat::FeaturePlot(SingleCellData,
                                  features = "Hal",
                                  pt.size = 1)+
    ggplot2::ggtitle("Hal - Portal Marker")+
    ggplot2::ylab("")+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 22))

#####Figure 4####

options(ggrepel.max.overlaps = Inf)

design_layout4 <- "
1122
3344
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
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20,
                                                    face = "plain"))

patchworktest4 <- patchworktest4 + patchwork::plot_annotation(title = "Figure 4",
                                            theme = theme(plot.title = ggplot2::element_text(size = 26)))
       grDevices::pdf(here::here("data/Figure4.pdf"), height = 20, width = 20)
        patchworktest4
        dev.off()

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
                            "CS2" = "CS",
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
SeuratHepatocytes$hash.mcl.ID <- Idents(SeuratHepatocytes)
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

#####Supporting figure 4####
#included from a reviewer suggestion

Seurat::Idents(SingleCellData)<- "hash.mcl.ID"

DimPlotIndividual_L1 <- Seurat::DimPlot(SingleCellData,
                                     pt.size = 0.7,
                                     label = T,
                                     label.size = 6,
                                     repel = T,
                                     # cols =  wesanderson::wes_palette(24,
                                     #                                 name = "Zissou1",
                                     #                                 type = "continuous"),
                                     cols = viridis::turbo(24,begin = 0.2, end  =0.8),
                                     label.box = T,
                                     split.by = "Group") +
    ggplot2::ggtitle("DimPlot By individual") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none",
        strip.text = ggplot2::element_text(size = 22)


    )
DimPlotIndividual_L1

patchwork_sup2 <- DimPlotIndividual_L1+ patchwork::plot_annotation(title = "Supporting Fig. 2",
                                                              theme = theme(plot.title = ggplot2::element_text(size = 26)))

grDevices::pdf(here::here("data/SupportingFig4.pdf"), height = 20, width = 20)
patchwork_sup2
dev.off()

#####SupportingFig5####
#this used to be figure 5. Hence annotation for this and subsequent figures
#is slightly mismatched. I was too lazy to fix it, and too scared of
#breaking dependencies to move things around. Please accept my apologies.

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
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20,
                                                    face = "plain"))

patchworktest5 <- patchworktest5 + patchwork::plot_annotation(title = "Supporting Fig. 3",
                                                                               theme = theme(plot.title = ggplot2::element_text(size = 26)))

    grDevices::pdf(here::here("data/SupportingFig5.pdf"), height = 20, width = 20)
     patchworktest5
     dev.off()

#####PseudoDEG on hepatocytes####
PseudoHepatocytes <- Pseudobulk(SeuratHepatocytes)
DEGHepatocytes <-DGEPseudo(PseudoHepatocytes)

#write excel file
write_excel_file(DEGHepatocytes, "SupportingFile8")
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

GOCC_ExcelExport <- list("PHvsLUp"=GOCCPHvsLUp@result,
                         "LvsPHUp"= GOCCLvsPHUp@result)
#write the merged object as Excel file
write_excel_file(GOCC_ExcelExport, "SupportingFile9")

#####Figure 5####
setup_figure_5  <- png::readPNG(here::here("data-raw/Figure 5 setup.png"))
setup_figure_5 <- grid::rasterGrob(grDevices::as.raster(setup_figure_5),interpolate = T)

design_layout6 <- "
1222
3344
3344
3344
"


patchworktest6 <- patchwork::free(
    patchwork::wrap_elements(setup_figure_5))+
    UpsetHepatocytes+
    GOCCLvsPHUpTreePlot+
    GOCCPHvsLUpTreePlot+
    patchwork::plot_layout(design = design_layout6)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))
patchworktest6 <- patchworktest6 + patchwork::plot_annotation(title = "Figure 5",
                                                                            theme = theme(plot.title = ggplot2::element_text(size = 26)))
   grDevices::pdf(here::here("data/Figure5.pdf"), height = 20, width = 20)
   patchworktest6
   dev.off()

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
    ggplot2::ggtitle("DimPlot by cluster \n for PH samples") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5,
                                           size = 24),
        legend.position = "none"
    )

#assign cell identities based on marker gene expression
Seurat::Idents(SeuratPH)<- "seurat_clusters"
SeuratPH <- RenameIdents(SeuratPH,
                            "0" = "Dedifferentiated \nHepatocytes",
                            "1" = "Dedifferentiated \nHepatocytes",
                            "2" = "Dedifferentiated \nHepatocytes",
                            "3" = "Dedifferentiated \nHepatocytes",
                            "4" = "Dedifferentiated \nHepatocytes",
                            "5" = "Dedifferentiated \nHepatocytes",
                            "6" = "Dedifferentiated \nEndothelial Cells",
                            "7" = "Dedifferentiated \nHepatocytes",
                            "8" = "Dedifferentiated \nHepatocytes",
                            "9" = "Dedifferentiated \nHepatocytes",
                            "10" = "Dedifferentiated \nHepatocytes",
                            "11" = "Dedifferentiated \nHepatocytes",
                            "12" = "Dedifferentiated \nHepatocytes",
                            "13" = "Dedifferentiated \nHepatocytes",
                            "14" = "Dedifferentiated \nHepatocytes",
                            "15" = "Dedifferentiated \nHepatocytes",
                            "16" = "Dedifferentiated \nHepatocytes",
                            "17" = "Stellate Cells",
                            "18" = "Dedifferentiated \nHepatocytes",
                            "19" = "Dedifferentiated \nHepatocytes",
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
                                                       size = 14),
                   axis.text.y = ggplot2::element_text(size = 18),
                   axis.title = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 24,
                                                      hjust = 0.5,
                                                      face = "bold"))

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
    ggplot2::ggtitle("Cell distribution in sample types") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(size = 24,
                                      hjust = 0.5,
                                      face = "bold"),
        axis.text.y = ggplot2::element_text(size = 16),
        axis.text.x = ggplot2::element_text(size = 14),
        axis.title.x = ggplot2::element_text(size = 16),
        legend.text = ggplot2::element_text(size = 12),
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "black")

    )+
    ggplot2::ylab("")+
    ggplot2::xlab("No. of cells")+
    ggplot2::scale_fill_manual(
        values = wesanderson::wes_palette("Darjeeling1",
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
    ggplot2::ggtitle("% cells in sample types") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(size = 24,
                                      hjust = 0.5,
                                      face = "bold"),
        axis.text.y = ggplot2::element_text(size = 16),
        axis.text.x = ggplot2::element_text(size = 14),
        axis.title.x = ggplot2::element_text(size = 16),
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


#####Figure 6####
design_layout7 <- "
11#222
11#222
11#222
333444
333444
333555
###555
667788
667788
"

setup_figure_6 <- png::readPNG(here::here("data-raw/Figure 6 setup.png"))
setup_figure_6 <- grid::rasterGrob(grDevices::as.raster(setup_figure_6),interpolate = T)

setup_figure_6B <- png::readPNG(here::here("data-raw/Figure 6B setup.png"))
setup_figure_6B<- grid::rasterGrob(grDevices::as.raster(setup_figure_6B),interpolate = T)

patchworktest7 <-
    patchwork::free(
        patchwork::wrap_elements(setup_figure_6)
    )+
    patchwork::free(DimplotCelltypePH, type = "label")+
    patchwork::free(DotplotFigurePH)+
    cellcountfigure+
    cellratio+
    patchwork::free(
        patchwork::wrap_elements(setup_figure_6B)
    )+
   patchwork::free(DimplotEndothelial, type = "label")+
    DotplotFigureEndo+
    patchwork::plot_layout(design = design_layout7)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20, face = "plain"))
patchworktest7 <- patchworktest7 + patchwork::plot_annotation(title = "Figure 6",
                                                              theme = theme(plot.title = ggplot2::element_text(size = 26,
                                                                                                               face = "plain")))
       grDevices::pdf(here::here("data/Figure6.pdf"), height = 20, width = 20)
       patchworktest7
       dev.off()

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

trimmed_cpm <- proteomicsDataAnno|>
    dplyr::select(-ACCID) |>
    tidyr::pivot_longer(cols = 1:24, values_to = "logAbundance",names_to = "SampleID")

#Select Candidates in Gene List


#create annotation key for heatmap
setup <- targets::tar_read(ProteomicsSetup)
key <- as.data.frame(setup)
key <- key |>
    dplyr::select(Tissue, SampleID)
key$Tissue <- factor(key$Tissue, c("liver", "CS", "PH"))

trimmed_cpm <- dplyr::filter(trimmed_cpm,
                             SYMBOL %in% candidates$gene)
trimmed_cpm <- left_join(trimmed_cpm, key)
trimmed_cpm$SampleID <- factor(trimmed_cpm$SampleID, levels = key$SampleID)

heatmap_candidates <- tidyHeatmap::heatmap(trimmed_cpm,
                                       .row = SYMBOL,
                                       .column = SampleID,
                                       .value = logAbundance,
                                       scale = "column",
                                       column_dend_height = unit(0, "cm"),
                                       cluster_rows = T,
                                       cluster_columns = F,
                                       palette_value = circlize::colorRamp2(c(-3,-2,-1,0,1,2,3), viridis::inferno(7)),
                                       row_names_gp = ggfun::gpar(fontsize = 16),
                                       column_title = "Candidates from proteomics data",
                                       column_title_gp = ggfun::gpar(fontsize = 24),
                                       show_row_names = T,
                                       show_column_names = F,
                                       show_heatmap_legend = T,
                                       row_title = NULL
) |> tidyHeatmap::annotation_tile(Tissue,palette = c("red", "cyan", "orange"),
                                  show_legend = T,
                                  show_title = FALSE,
                                  annotation_name = NULL)

#extract heatmap order
as_complex_heatmap <- tidyHeatmap::as_ComplexHeatmap(heatmap_candidates)
as_complex_heatmap <- ComplexHeatmap::draw(as_complex_heatmap)
row_order_heatmap <- ComplexHeatmap::row_order(as_complex_heatmap)
row_order_symbols <- trimmed_cpm |>
    dplyr::distinct(SYMBOL) |>
    dplyr::arrange(SYMBOL) |>
    dplyr::pull()
row_order_symbols <- row_order_symbols[row_order_heatmap]

#repeat group UMAP for identification
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


glut1Plot <- Seurat::FeaturePlot(SingleCellData,
                         features = "Slc2a1",
                         pt.size = 1)+
    ggplot2::ggtitle("Glut1 Expression")+
    ggplot2::theme(title = ggplot2::element_text(size = 22))

glut2Plot <- Seurat::FeaturePlot(SingleCellData,
                                 features = "Slc2a2",
                                 pt.size = 1)+
    ggplot2::ggtitle("Glut2 Expression")+
    ggplot2::theme(title = ggplot2::element_text(size = 22))

Idents(SingleCellData)<-"celltype"
DotplotCandidateFigure <- Seurat::DotPlot(SingleCellData,
                                 features = rev(row_order_symbols))+
    ggplot2::ggtitle("Candidate Genes - RNA expression")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       vjust = 0.65,
                                                       size = 18),
                   axis.text.y = ggplot2::element_text(size = 18),
                   plot.title = ggplot2::element_text(hjust = 0.5,
                                                      size = 24)
                   )+
    ggplot2::coord_flip()


#####Figure 7####
design_layout8 <- "
111222
111222
334455
"
patchworktest8 <-
    tidyHeatmap::wrap_heatmap(heatmap_candidates)+
    DotplotCandidateFigure+
    patchwork::free(DimPlotGroup)+
    patchwork::free(glut1Plot)+
    patchwork::free(glut2Plot)+
    patchwork::plot_layout(design = design_layout8)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20,
                                                    face = "plain"))

patchworktest8 <- patchworktest8 + patchwork::plot_annotation(title = "Figure 7",
                                                              theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 26)))
      grDevices::pdf(here::here("data/Figure7.pdf"), height = 20, width = 20)
      patchworktest8
      dev.off()





