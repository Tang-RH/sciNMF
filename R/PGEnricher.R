#' Program Enricher
#'
#' Performing universal enrichment analysis on transcriptional programs (gene sets)
#'
#' @param pg.list A list containing transcriptional programs, each program is composed with genes.
#' @param gene.set User input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene. All built-in data.frame object
#' are \code{"GO_BP", "GO_CC", "GO_MM", "Hallmark", "CellMarker"}. You can set a character vector for using the built-in gene sets. Default is "GO_BP"
#' @param universe background genes. If missing, the all genes listed in the database (eg gene.set data.frame) will be used as background.
#' @param only.df Logical, if TRUE, the function will return only the data.frame for each enrichment results.
#' If FALSE, the full enrichment result object will be returned. Default is TRUE.
#' @param minGSSize The minimum number of genes required for a gene set to be considered in the analysis.
#' Default is 10.
#' @param maxGSSize The maximum number of genes allowed in a gene set to be considered in the analysis.
#' Default is 500.
#' @param padj.cutoff The adjusted p-value cutoff for selecting significantly enrichment results.
#' Default is 0.05.
#' @param padj.method The method used for multiple testing correction of p-values.
#' Options include "holm", "hochberg", "BH", and other methods supported by \code{\link[stats]{p.adjust}}.
#' Default is "BH" (Benjamini-Hochberg).
#' @param ... Additional arguments to be passed to the \code{\link[clusterProfiler]{enricher}} of clusterProfiler package.
#'
#' @return A list containing the enrichment results for each input transcriptional program.
#' Each element of the list represents one transcriptional program and contains a data.frame with enrichment results.
#'
#' @details
#' The function uses the \code{clusterProfiler::\link[clusterProfiler]{enricher}} function from the clusterProfiler package to perform gene set enrichment analysis based on \code{stats::\link[stats]{phyper}} hypergeometric test.
#' It calculates the enrichment of programs in the provided gene list using the specified gene set database.
#' The resulting data frame contains information about enriched gene sets, including gene set ID, name, p-value,
#' adjusted p-value, and other statistics.
#'
#'
#' @seealso
#' \code{\link[clusterProfiler]{enricher}} in the clusterProfiler package for more details of the enrichment analysis.
#'
#' \code{\link[gson]{read.gmt}} for gene.set data.frame generation.
#'
#' \code{\link{ls_gs_sciNMF}} built-in gene sets.
#'
#' \url{https://maayanlab.cloud/Enrichr/#libraries} The source of gene.set.
#'
#' @importFrom clusterProfiler enricher
#' @importFrom utils  data
#'
#' @export
#'
PGEnricher = function(pg.list, gene.set = "GO_BP", universe = NULL,
                    only.df = TRUE,
                    minGSSize = 10, maxGSSize = 500,
                    padj.cutoff =0.05, padj.method = 'BH',...){
    # data("ls_gs_sciNMF", package = "sciNMF")

    if(is.character(gene.set)){
        if(!all(gene.set %in% names(ls_gs_sciNMF))){
            warning(paste(setdiff(gene.set, names(ls_gs_sciNMF)), collapse = ', '), '  not supported')
            gene.set = intersect(gene.set, names(ls_gs_sciNMF))
            if(length(gene.set)<1){stop('Please set correct gene.set')}
        }
        TERM2GENE = do.call(what = rbind, c(ls_gs_sciNMF[gene.set],make.row.names = FALSE))
    }else{
        TERM2GENE = gene.set
    }

    ls_res_enrich = lapply(pg.list,function(pg){
        res = clusterProfiler::enricher(pg, TERM2GENE = TERM2GENE,
                                        minGSSize = minGSSize,
                                        maxGSSize = maxGSSize,
                                        pAdjustMethod = padj.method,
                                        universe = universe,
                                        pvalueCutoff = padj.cutoff,...)
        if(only.df){
            res = res@result
        }
        return(res)
    })
    return(ls_res_enrich)
}
