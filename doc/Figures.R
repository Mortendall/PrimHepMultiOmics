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

setup_figure_ggplot <- ggplot2::ggplot()+ggplot2::annotation_custom(setup_figure,
                                                                    xmin = -Inf,
                                                                    xmax = Inf,
                                                                    ymin = -Inf,
                                                                    ymax = Inf)


#decided to remove cpm-figs from main text as they were mainly for internal QC

design_layout1 <- "
112333
112333
444555
444555
667788
667788
66##88
66####
"
patchworktest <- patchwork::wrap_elements(setup_figure)+
    pcaplotRNA +
    upsetRNA +
    #cpm_figs +
    upsetGO+
    CompareClusterFigure_All+
    MitoHeatmapRNA+
    RiboHeatmapRNA+
    ECMHeatmapRNA+
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

setup_figure_ggplot <- ggplot2::ggplot()+ggplot2::annotation_custom(setup_figure,
                                                                    xmin = -Inf,
                                                                    xmax = Inf,
                                                                    ymin = -Inf,
                                                                    ymax = Inf)

#decided to remove cpm-figs from main text as they were mainly for internal QC

design_layout2 <- "
112333
112333
444555
444555
667788
667788
66##88
66####
"
patchworktest2 <- patchwork::wrap_elements(setup_figure_2)+
    pcaplotProt +
    upsetProt +
    #cpm_figsProt +
    upsetGOProt+
    CompareClusterFigure_All_Prot+
    HeatmapProteomeMito+
    HeatmapProteomeRibo+
    HeatmapProteomeECM+
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

setup_figure_ggplot <- ggplot2::ggplot()+ggplot2::annotation_custom(setup_figure,
                                                                    xmin = -Inf,
                                                                    xmax = Inf,
                                                                    ymin = -Inf,
                                                                    ymax = Inf)

design_layout3 <- "
1233
4455
4455
6655
6655
"

patchworktest3 <- patchwork::wrap_elements(setup_figure_3)+
    (mdsPseudo+ggplot2::theme(axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = -200, unit = "pt"))))+
    UpsetPseudo+
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
                            here::here("doc/hepamorphosis_fig.pdf")),
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
