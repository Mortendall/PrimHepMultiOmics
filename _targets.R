# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tidyverse)
library(Seurat)
library(SummarizedExperiment)
library(org.Mm.eg.db)
# library(tarchetypes) # Load other packages as needed. # nolint
# Set target options:
tar_option_set(
  packages = desc::desc_get_deps()$package[-1],
  format = "rds"
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options. # nolint

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
list(
  tar_target(
    name = bulkLoad,
    command = count_matrix_assembly("count_matrix.xlsx")
  ),
  tar_target(
    name = metadataLoad,
    command = load_metadata("metadata.xlsx")
  ),
  tar_target(
    name = filteredMetadata,
    command = removeOutlierMeta(metadataLoad, c("558L", "544CS", "544PH"))
  ),
  tar_target(
    name = filteredCounts,
    command = removeOutlierCounts(bulkLoad, filteredMetadata)
  ),
  tar_target(
    name = designMatrix,
    command = Generate_design_matrix(filteredMetadata)
  ),
  tar_target(
    name = ctrsts,
    command = Generate_ctrst(designMatrix)
  ),
  tar_target(
    name = DGElist,
    command = DGEprep(filteredMetadata, filteredCounts, designMatrix)
  ),
  tar_target(
    name = DGETable,
    command = DGEanalysis(DGElist, designMatrix, ctrsts)
  ),
  tar_target(
    name = MDSplot,
    command = MDSanalysis(DGElist, filteredMetadata)
  ),
  tar_target(
    name = upsetplot,
    command = UpsetplotGeneration(DGETable)
  ),
  tar_target(
    name = cpmPlot,
    command = XYCPMplot(DGElist)
  ),
  tar_target(
    name = GOCC,
    command = GOCCSplit(DGETable)
  ),
  tar_target(
    name = CCDotplot,
    command = DotplotCC(GOCC,
                        "GSE analysis - Upregulated and downregulated\n RNAseq analysis")
  ),
  tar_target(
    name = UpsetplotGO,
    command = UpsetGO(GOCC)
  ),
  tar_target(
      name = HeatmapMitoRNA,
      command= RNAHeatmap(GOCC,43,DGElist,filteredMetadata, T)
  ),
  tar_target(
      name = WriteDEGTable,
      command = write_excel_file(DGETable, "SupportingFile1")
  ),
  tar_target(
      name = WriteGOtable,
      command = write_excel_file(GOCC@compareClusterResult, "SupportingFile2")
  ),
  tar_target(
    name = RawProteomics,
    command = ProteomicsDataLoader("/morten_proteome_combined_hybrid.txt")
  ),
  tar_target(
    name = ConversionKey,
    command = IDkey(RawProteomics)
  ),
  tar_target(
    name = ProteomicsSetup,
    command = GenerateSetup(RawProteomics)
  ),
  tar_target(
    name = NormalizedMatrix,
    command = DataFiltering(
      RawProteomics,
      ProteomicsSetup
    )
  ),
  tar_target(
    name = DAPResults,
    command = LimmaAnalysis(
      NormalizedMatrix,
      ProteomicsSetup,
      ConversionKey
    )
  ),
  tar_target(
      name = DAPRexcel,
      command = write_excel_file(DAPResults, "SupportingFile3")
  ),
  tar_target(
    name = MDSproteomics,
    command = ProteomicsMDS(
      NormalizedMatrix,
      ProteomicsSetup
    )
  ),
  tar_target(
    name = GOCCProteomics,
    command = GOCCSplit(DAPResults)
  ),
  tar_target(
      name = GOCCExcelProt,
      command = write_excel_file(GOCCProteomics@compareClusterResult,
                                 "SupportingFile4")
  ),
  tar_target(
    name = DotplotCCProteomics,
    command = DotplotCC(GOCCProteomics,
                        "GSE analysis - Upregulated and downregulated")
  ),
  tar_target(
    name = UpsetplotGOProteomics,
    command = UpsetGO(GOCCProteomics)
  ),
  tar_target(
    name = CPMplotProteomics,
    command = CPMPlotsProteomics(NormalizedMatrix)
  ),
  tar_target(
    name = UpsetplotProteomics,
    command = UpsetProteomics(DAPResults)
  ),
  tar_target(
    name = HeatmapProteome,
    command = proteomicsHeatmap(GOCCProteomics,
                                targetrow = 100,
                                counts = NormalizedMatrix,
                                ProteomicsSetup,
                                T)
  ),
  tar_target(
    name = seuratObject,
    command = loadSingleCell("220503_liver_full-seurat_updated.rds")
  ),
  tar_target(
    name = pseudoDEseq,
    command = Pseudobulk(seuratObject)
  ),
  tar_target(
    name = PseudoCounts,
    command = DGECountMatrix(pseudoDEseq)
  ),
  tar_target(
    name = PseudoDEG,
    command = DGEPseudo(pseudoDEseq)
  ),
  tar_target(
      name = PseudoDEGExcel,
      command = write_excel_file(PseudoDEG,
                                 "SupportingFile5")
  ),
  tar_target(
    name = ProtRNACor,
    command = ProteinRNACorrelation(
      PseudoDEG,
      DAPResults)
  ),
  tar_target(
      name = CorrelationGOCC,
      command = GOCCRNAProt(ProtRNACor,
                            DAPResults)
  ),
  tar_target(
      name = CorrelationGOCCExcel,
      command = write_excel_file(CorrelationGOCC@compareClusterResult,
                                 "SupportingFile7")
  ),
  tar_target(
        name = CorrelationFigure,
        command = ProteinRNACorFigure(
            ProtRNACor
            )
  ),
  tar_target(
      name = CorrelationGOCCFigure,
      command = DotplotCC(CorrelationGOCC,
                          "GSE analysis - L vs PH")
  ),
  tar_target(
      name = PseudoMDS,
      command = MDSPseudoPlot(PseudoCounts)
  ),
  tar_target(
      name = GOCCPseudo,
      command = GOCCSplitPseudo(PseudoDEG, "CC")
  ),
  tar_target(
      name = GOCCPseudoExcel,
      write_excel_file(GOCCPseudo@compareClusterResult,
                       "SupportingFile6")
  ),
  tar_target(
      name = dotplotPseudo,
      command = DotplotCC(GOCCPseudo,
                          "GSE analysis - Upregulated and downregulated")
  ),
  tar_target(
      name = UpsetPseudo,
      command = UpsetplotGenerationPseudo(PseudoDEG,"Upset plot - Pseudobulk")
  ),
  #Add two more heatmaps to the pipeline
  tar_target(
      name = HeatmapECMRNA,
      command = RNAHeatmap(GOCC,56, DGElist, filteredMetadata, T)
  ),
  tar_target(
      name = HeatmapRiboRNA,
      command= RNAHeatmap(GOCC,102,DGElist,filteredMetadata, T)
  ),
  tar_target(
      name = HeatmapProteomeRibo,
      command = proteomicsHeatmap(GOCCProteomics,
                                  targetrow = 134,
                                  counts = NormalizedMatrix,
                                  ProteomicsSetup,
                                  T)
  ),
  tar_target(
      name = HeatmapProteomeECM,
      command = proteomicsHeatmap(GOCCProteomics,
                                  targetrow = 72,
                                  counts = NormalizedMatrix,
                                  ProteomicsSetup,
                                  T)
  ),
  tar_target(
      name = VolcanoRNA,
      command = volcano_plotter(DGETable, "RNA")
  ),
  tar_target(
      name = VolcanoProtein,
      command = volcano_plotter(DAPResults, "proteomics")
  ),
  tar_target(
      name = VolcanoPseudo,
      command = volcano_plotter(PseudoDEG, "pseudo")
  ),
  #make clustering
  tar_target(
      name = go_network_cc,
      command = prepare_network()
  ),
  #cluster mRNA data
  tar_target(
      name = mRNA_GO_data,
      command = prepare_data(GOCC)
  ),
  tar_target(
      name = mRNA_subgraphs,
      command = create_subgraphs(mRNA_GO_data, go_network_cc)
  ),
  tar_target(
      name = cluster_mRNA,
      command = cluster_go_terms(mRNA_subgraphs, go_network_cc, mRNA_GO_data, 0.6)
  ),
  tar_target(
      name = mRNA_plot_prepare,
      command = prepare_plot(cluster_mRNA,mRNA_subgraphs)
  ),
  tar_target(
      name = mRNA_plots,
      command = make_plots(cluster_mRNA, mRNA_plot_prepare)
  ),
  #cluster proteomics data
  tar_target(
      name = prot_GO_data,
      command = prepare_data(GOCCProteomics)
  ),
  tar_target(
      name = prot_subgraphs,
      command = create_subgraphs(prot_GO_data, go_network_cc)
  ),
  tar_target(
      name = cluster_prot,
      command = cluster_go_terms(prot_subgraphs, go_network_cc, prot_GO_data, 0.6)
  ),
  tar_target(
      name = prot_plot_prepare,
      command = prepare_plot(cluster_prot,prot_subgraphs)
  ),
  tar_target(
      name = prot_plots,
      command = make_plots(cluster_prot, prot_plot_prepare)
  )
)

