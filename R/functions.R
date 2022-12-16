#' count_matrix_loader
#'
#' @param file_type your raw data file type (presently optimized for txt)
#'
#' @return a count matrix containing raw count values for the experiment with group annotations.

count_matrix_assembly <- function(file_type){
    count_file <- fs::dir_ls(here::here("data-raw/"),
                             regexp = file_type,
                             recurse = TRUE)
    count_matrix <- openxlsx::read.xlsx(xlsxFile = count_file,sheet = 1)
    count_matrix$Geneid<-stringr::str_remove_all(count_matrix$Geneid, "\\..*")
    rownames(count_matrix) <- count_matrix$Geneid
    count_matrix <- count_matrix %>%
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
                            recurse = T)
    metadata <- openxlsx::read.xlsx(xlsxFile = data_file)
    metadata <- metadata |>
        dplyr::mutate(Group = dplyr::case_when(
            Group == "WT_L"~"L",
            Group == "WT_CS"~"CS",
            Group == "WT_PH"~"PH"
    return(metadata)
}
