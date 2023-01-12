#####Figure 1####
pcaplotRNA <- targets::tar_read(MDSplot)
upsetRNA <- targets::tar_read(upsetplot)
cpm_figs <- targets::tar_read(cpmPlot)
upsetGO <- targets::tar_read(UpsetplotGO)
CompareClusterFigure_All <- targets::tar_read(GOCCdotplot)

design_layout1 <- "
1222
3333
4455
4455
##55
"
patchworktest <- pcaplotRNA +
    upsetRNA +
    cpm_figs +
    upsetGO+
    CompareClusterFigure_All+
    patchwork::plot_layout(design = design_layout1)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

 # grDevices::pdf(here::here("data/Figure1.pdf"), height = 20, width = 20)
 # patchworktest
 # dev.off()

#####Figure 2####
pcaplotProt <- targets::tar_read(MDSproteomics)
upsetProt <- targets::tar_read(UpsetplotProteomics)
cpm_figsProt <- targets::tar_read(CPMplotProteomics)
upsetGOProt <- targets::tar_read(UpsetplotGOProteomics)
CompareClusterFigure_All_Prot <- targets::tar_read(DotplotCCProteomics)

design_layout2 <- "
1222
3333
4455
4455
##55
"
patchworktest2 <- pcaplotProt +
    upsetProt +
    cpm_figsProt +
    upsetGOProt+
    CompareClusterFigure_All_Prot+
    patchwork::plot_layout(design = design_layout2)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

# grDevices::pdf(here::here("data/Figure2.pdf"), height = 20, width = 20)
#  patchworktest2
#  dev.off()

#####Figure 3####
mdsPseudo <- targets::tar_read(PseudoMDS)
UpsetPseudo <- targets::tar_read(UpsetPseudo)
GOPseudo <- targets::tar_read(dotplotPseudo)
LogFCCompare <- targets::tar_read(CorrelationFigure)
GOcor <- targets::tar_read(CorrelationGOCCFigure)

design_layout3 <- "
1222
3344
3355
3355
##55
"

patchworktest3 <- mdsPseudo+
    UpsetPseudo+
    GOPseudo+
    LogFCCompare+
    GOcor+
    patchwork::plot_layout(design = design_layout3)+
    patchwork::plot_annotation(tag_levels = "A")&
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 20))

 # grDevices::pdf(here::here("data/Figure3.pdf"), height = 20, width = 20)
 # patchworktest3
 # dev.off()
