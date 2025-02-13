#' List of sciNMF built-in gene sets data.frame
#'
#' A list of annotation of TERM TO GENE mapping data.frames, which have 2 columns with term and gene. The data.frames are generated by \code{\link[gson]{read.gmt}}
#'
#' @usage ls_gs_sciNMF[["GO_BP"]]
#' @format
#' A list with 5 TERM TO GENE mapping data.frames and their sources:
#' \itemize{
#'   \item\code{CellMarker}: CellMarker_Augmented_2021.txt
#'   \item\code{GO_BP}: GO_Biological_Process_2023.txt
#'   \item\code{GO_CC}: GO_Cellular_Component_2023.txt
#'   \item\code{GO_MF}: GO_Molecular_Function_2023.txt
#'   \item\code{Hallmark}: MSigDB_Hallmark_2020.txt
#' }
#' @examples data("ls_gs_sciNMF")
#' @seealso \code{\link[gson]{read.gmt}} for gene.set data.frame generation.
#' @source \url{https://maayanlab.cloud/Enrichr/#libraries}
"ls_gs_sciNMF"



#' A seurat object with 400 cells and 13,642 features
#'
#' This object is used to test the package installation. It contains 400 leukemic cells from 4 AML patients in our study.
#'
#' @format
#' \describe{
#'   \item{Seurat object}{
#'     13642 features across 400 samples within 1 assay
#'     Active assay: RNA (13642 features, 0 variable features)
#'   }
#' }
#' @examples data("SrtObj")
"SrtObj"