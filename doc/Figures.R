#####Figure 1####
pcaplotRNA <- targets::tar_read(MDSplot)
upsetRNA <- targets::tar_read(upsetplot)
cpm_figs <- targets::tar_read(cpmPlot)
upsetGO <- targets::tar_read(UpsetplotGO)
CompareClusterFigure_All <- targets::tar_read(CCDotplot)
MitoHeatmapRNA <- targets::tar_read(HeatmapMitoRNA)
RiboHeatmapRNA <- targets::tar_read(HeatmapRiboRNA)
ECMHeatmapRNA <- targets::tar_read(HeatmapECMRNA)
setup_figure <- png::readPNG(here::here("data-raw/Figure 1 setup.png"))
setup_figure <- grid::rasterGrob(grDevices::as.raster(setup_figure),interpolate = T)

# setup_figure_ggplot <- ggplot2::ggplot()+ggplot2::annotation_custom(setup_figure,
#                                                                     xmin = -Inf,
#                                                                     xmax = Inf,
#                                                                     ymin = -Inf,
#                                                                     ymax = Inf)


#decided to remove cpm-figs from main text as they were mainly for internal QC

# design_layout1 <- "
# 112333
# 112333
# 444555
# 444555
# 667788
# 667788
# 66##88
# 66####
# "
 design_layout1 <- "
 112333
 112333
 444555
 444555
 666666
 777777
 888888
 ######
 "

patchworktest <- patchwork::free(patchwork::wrap_elements(setup_figure))+
    pcaplotRNA +
    upsetRNA +
    #cpm_figs +
    upsetGO+
    CompareClusterFigure_All+
    tidyHeatmap::wrap_heatmap(MitoHeatmapRNA)+
    tidyHeatmap::wrap_heatmap(RiboHeatmapRNA)+
    tidyHeatmap::wrap_heatmap(ECMHeatmapRNA)+
    patchwork::plot_layout(design = design_layout1)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

patchworktest <- patchworktest + patchwork::plot_annotation(title = "Figure 1",
                                                            theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 26)))
      grDevices::pdf(here::here("data/Figure1.pdf"), height = 20, width = 20)
      patchworktest
      dev.off()

#####Figure 2####
pcaplotProt <- targets::tar_read(MDSproteomics)
upsetProt <- targets::tar_read(UpsetplotProteomics)
cpm_figsProt <- targets::tar_read(CPMplotProteomics)
upsetGOProt <- targets::tar_read(UpsetplotGOProteomics)
CompareClusterFigure_All_Prot <- targets::tar_read(DotplotCCProteomics)
HeatmapProteomeMito <- targets::tar_read(HeatmapProteome)
HeatmapProteomeRibo <- targets::tar_read(HeatmapProteomeRibo)
HeatmapProteomeECM <- targets::tar_read(HeatmapProteomeECM)
setup_figure_2 <- png::readPNG(here::here("data-raw/Figure 2 setup.png"))
setup_figure_2 <- grid::rasterGrob(grDevices::as.raster(setup_figure_2),interpolate = T)

# setup_figure_ggplot_2 <- ggplot2::ggplot()+ggplot2::annotation_custom(setup_figure_2,
#                                                                     xmin = -Inf,
#                                                                     xmax = Inf,
#                                                                     ymin = -Inf,
#                                                                     ymax = Inf)

#decided to remove cpm-figs from main text as they were mainly for internal QC

design_layout2 <- "
 112333
 112333
 444555
 444555
 666666
 777777
 888888
 ######
 "
patchworktest2 <- patchwork::free(patchwork::wrap_elements(setup_figure_2))+
    pcaplotProt +
    upsetProt +
    #cpm_figsProt +
    upsetGOProt+
    CompareClusterFigure_All_Prot+
    tidyHeatmap::wrap_heatmap(HeatmapProteomeMito)+
    tidyHeatmap::wrap_heatmap(HeatmapProteomeRibo)+
    tidyHeatmap::wrap_heatmap(HeatmapProteomeECM)+
    patchwork::plot_layout(design = design_layout2)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

patchworktest2 <- patchworktest2 + patchwork::plot_annotation(title = "Figure 2",
                                                                theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 26)))


      grDevices::pdf(here::here("data/Figure2.pdf"), height = 20, width = 20)
       patchworktest2
       dev.off()

#####Figure 3####
mdsPseudo <- targets::tar_read(PseudoMDS)
UpsetPseudo <- targets::tar_read(UpsetPseudo)
GOPseudo <- targets::tar_read(dotplotPseudo)
LogFCCompare <- targets::tar_read(CorrelationFigure)
GOcor <- targets::tar_read(CorrelationGOCCFigure)
setup_figure_3 <- png::readPNG(here::here("data-raw/Figure 3 setup.png"))
setup_figure_3 <- grid::rasterGrob(grDevices::as.raster(setup_figure_3),interpolate = T)



design_layout3 <- "
112233
112233
444555
444555
666555
666555
"

patchworktest3 <- patchwork::free(patchwork::wrap_elements(setup_figure_3))+
    (mdsPseudo+ggplot2::theme(
     #   axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = -200, unit = "pt"))
        )
     )+
    patchwork::free(UpsetPseudo)+
    GOPseudo+
    (LogFCCompare+ggplot2::theme(axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 0, unit = "pt")),
                                 axis.title.x = ggplot2::element_text(vjust = 5)))+
    GOcor+
    patchwork::plot_layout(design = design_layout3)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

patchworktest3 <- patchworktest3 + patchwork::plot_annotation(title = "Figure 3",
                                                              theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 26)))

    grDevices::pdf(here::here("data/Figure3.pdf"), height = 20, width = 20)
    patchworktest3
    dev.off()

    #####Supporting fig 1####
    RNA_volcano <- targets::tar_read(VolcanoRNA)
    RNA_clustering <- targets::tar_read(mRNA_plots)
    design_layout_sup1 <- "
112233
444555
666777
888999
"

    patchwork_sup1 <- patchwork::free(RNA_volcano[[1]])+
        patchwork::free(RNA_volcano[[2]])+
        patchwork::free(RNA_volcano[[3]])+
        patchwork::free(RNA_clustering[[1]])+
        patchwork::free(RNA_clustering[[2]])+
        patchwork::free(RNA_clustering[[3]])+
        patchwork::free(RNA_clustering[[4]])+
        patchwork::free(RNA_clustering[[5]])+
        patchwork::free(RNA_clustering[[6]])+
        patchwork::plot_layout(design = design_layout_sup1)+
        patchwork::plot_annotation(tag_levels = "A")&
        ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

    patchwork_sup1 <- patchwork_sup1 + patchwork::plot_annotation(title = "Supporting Fig. 1",
                                                                  theme = theme(plot.title = ggplot2::element_text(size = 26)))

    grDevices::pdf(here::here("data/SupportingFig1.pdf"), height = 20, width = 20)
    patchwork_sup1
    dev.off()

    #####Supporting Figure 2####

    Protein_volcano <- targets::tar_read(VolcanoProtein)
    Protein_clustering <- targets::tar_read(prot_plots)
    design_layout_sup2 <- "
112233
444555
666777
888999
"

    patchwork_sup2 <- patchwork::free(Protein_volcano[[1]])+
        patchwork::free(Protein_volcano[[2]])+
        patchwork::free(Protein_volcano[[3]])+
        patchwork::free(Protein_clustering[[1]])+
        patchwork::free(Protein_clustering[[2]])+
        patchwork::free(Protein_clustering[[3]])+
        patchwork::free(Protein_clustering[[4]])+
        patchwork::free(Protein_clustering[[5]])+
        patchwork::free(Protein_clustering[[6]])+
        patchwork::plot_layout(design = design_layout_sup1)+
        patchwork::plot_annotation(tag_levels = "A")&
        ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

    patchwork_sup2 <- patchwork_sup2 + patchwork::plot_annotation(title = "Supporting Fig. 2",
                                                                  theme = theme(plot.title = ggplot2::element_text(size = 26)))

    grDevices::pdf(here::here("data/SupportingFig2.pdf"), height = 20, width = 20)
    patchwork_sup2
    dev.off()


    #####Supporting Figure 3####
    #Used to be Figure 9 is a GO analysis of the terms that do not change in hepatocytes.
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
        ggplot2::ggtitle("GSE Analysis for Molecular function \nProteins with unchanged abundance between L and PH")+
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 24,
                                               hjust = 0.5))
    NormalizedMatrix <- targets::tar_read(NormalizedMatrix)
    ProteomicsSetup <- targets::tar_read(ProteomicsSetup)
    #prepare 3 heatmaps for figure. One for catalytic activity (term 2), one for
    # glutatione transferase activity (term 1) and one for peptidase activity (term 7)
    CatActHm <- proteomicsHeatmap(unchangedProtMF,
                                  targetrow = 2,
                                  counts = NormalizedMatrix,
                                  ProteomicsSetup,
                                  T)
    glutathioneHm <- proteomicsHeatmap(unchangedProtMF,
                                       targetrow = 1,
                                       counts = NormalizedMatrix,
                                       ProteomicsSetup,
                                       T)
    peptidaseHm <- proteomicsHeatmap(unchangedProtMF,
                                     targetrow = 7,
                                     counts = NormalizedMatrix,
                                     ProteomicsSetup,
                                     T)

    design_layout9 <- "
1122
1122
33##
44##
"

    patchworktest9 <-UnchangedDotPlot+
        tidyHeatmap::wrap_heatmap(CatActHm)+
        tidyHeatmap::wrap_heatmap(glutathioneHm)+
        tidyHeatmap::wrap_heatmap(peptidaseHm)+
        patchwork::plot_layout(design = design_layout9)+
        patchwork::plot_annotation(tag_levels = "A")&
        ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

    patchworktest9 <- patchworktest9 + patchwork::plot_annotation(title = "Supporting Fig. 2",
                                                                  theme = theme(plot.title = ggplot2::element_text(size = 26)))

    grDevices::pdf(here::here("data/SupportingFig3.pdf"), height = 20, width = 20)
    patchworktest9
    dev.off()


#####Supplemental figure 6####
design_supplemental4 <- "
    12
    3#
    "

    hepafigure_1 <- png::readPNG(here::here("data-raw/hepamorphosis_1.png"))
    hepafigure_1 <- grid::rasterGrob(grDevices::as.raster(hepafigure_1),
                                     interpolate = T)
    hepafigure_2 <- png::readPNG(here::here("data-raw/hepamorphosis_2.png"))
    hepafigure_2 <- grid::rasterGrob(grDevices::as.raster(hepafigure_2),
                                     interpolate = T)
    hepafigure_3 <- png::readPNG(here::here("data-raw/hepamorphosis_3.png"))
    hepafigure_3 <- grid::rasterGrob(grDevices::as.raster(hepafigure_3),
                                     interpolate = T)
    patchwork_figure_sup4 <- patchwork::free(patchwork::wrap_elements(hepafigure_1))+
        patchwork::free(patchwork::wrap_elements(hepafigure_2))+
        patchwork::free(patchwork::wrap_elements(hepafigure_3))+
        patchwork::plot_layout(design = design_supplemental4)+
        patchwork::plot_annotation(tag_levels = "A")&
        ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

    patchwork_figure_sup4 <- patchwork_figure_sup4+
        patchwork::plot_annotation(title = "Supporting Figure 4",
                                   theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 26)))

    grDevices::pdf(here::here("data/SupportingFig6.pdf"), height = 20, width = 20)
    patchwork_figure_sup4
    dev.off()



    #####Figure 4-8####

#as target is not well-suited to handle Seurat objects due to size,
#figure 4-9 are generated in separate script. Code is available in
#here::here("doc/SingleCellFigures.R")

#####Patch figures together####
qpdf::pdf_combine(input = c(here::here("data/Figure1.pdf"),
                            here::here("data/Figure2.pdf"),
                            here::here("data/Figure3.pdf"),
                            here::here("data/Figure4.pdf"),
                            here::here("data/Figure5.pdf"),
                            here::here("data/Figure6.pdf"),
                            here::here("data/Figure7.pdf"),
                            here::here("data/Figure8.pdf"),
                            here::here("data/SupportingFig1.pdf"),
                            here::here("doc/hepamorphosis_fig.pdf")),
                  output = here::here("data/Figures.pdf"))

#combined file is too large. Make two files instead

qpdf::pdf_combine(input = c(here::here("data/Figure5.pdf"),
                            here::here("data/Figure6.pdf"),
                            here::here("data/Figure7.pdf"),
                            here::here("data/SupportingFig1.pdf"),
                            here::here("data/SupportingFig2.pdf"),
                            here::here("data/SupportingFig3.pdf"),
                            here::here("data/SupportingFig4.pdf"),
                            here::here("data/SupportingFig5.pdf"),
                            here::here("data/SupportingFig6.pdf")),
                  output = here::here("data/Figures5_sup.pdf"))
qpdf::pdf_combine(input = c(here::here("data/Figure1.pdf"),
                            here::here("data/Figure2.pdf"),
                            here::here("data/Figure3.pdf"),
                            here::here("data/Figure4.pdf")),
                  output = here::here("data/Figures1_4.pdf"))

#calculate how many proteins were identified pr. sample type

proteomicsdata<- targets::tar_read(RawProteomics)
liver <- proteomicsdata[1:8] |>
    dplyr::rowwise() |>
    dplyr::summarise(rowsum = sum(dplyr::c_across(dplyr::everything()), na.rm = T)) |>
    dplyr::filter(rowsum!=0)
cs <- proteomicsdata[9:16]|>
    dplyr::rowwise() |>
    dplyr::summarise(rowsum = sum(dplyr::c_across(dplyr::everything()), na.rm = T)) |>
    dplyr::filter(rowsum!=0)
ph<- proteomicsdata[17:24]|>
    dplyr::rowwise() |>
    dplyr::summarise(rowsum = sum(dplyr::c_across(dplyr::everything()), na.rm = T)) |>
    dplyr::filter(rowsum!=0)
