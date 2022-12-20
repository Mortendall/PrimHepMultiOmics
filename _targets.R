# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = desc::desc_get_deps()$package[-1],
  format = "rds"

)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

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
    command = removeOutlierMeta(metadataLoad, "558L")
  ),
  tar_target(
    name = filteredCounts,
    command = removeOutlierCounts(bulkLoad, filteredMetadata, "558L")
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
      name = DGETable,
      command = DGEanalysis(filteredCounts, filteredMetadata, designMatrix, ctrsts)
  )
)
