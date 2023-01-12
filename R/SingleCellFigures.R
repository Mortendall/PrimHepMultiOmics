#The UMAP plots are too large to save as target objects. Instead, they are
#generated and saved through this script

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
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       vjust = 0.8,
                                                       size = 18),
                   axis.text.y = ggplot2::element_text(size = 18))


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

#####Make Cyp2e1 and Cyp2f2 plots####
Cyp2e1Plot <- Seurat::FeaturePlot(SingleCellData,
                                  features = "Cyp2e1",
                                  pt.size = 1)+
    ggplot2::ggtitle("Cyp2e1 - Central Marker")+
    ggplot2::theme(title = ggplot2::element_text(size = 22))

Cyp2f2Plot <- Seurat::FeaturePlot(SingleCellData,
                                  features = "Cyp2f2",
                                  pt.size = 1)+
    ggplot2::ggtitle("Cyp2f2 - Portal Marker")+
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


