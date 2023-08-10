#' Construct Meta Programs from NMF and Clustering Results
#'
#' This function extracts meta-programs from the results of Non-negative Matrix Factorization (NMF) and cluster relationship.
#'
#' @param WH.list
#' A nested list containing NMF results for each individual, generated from the \code{\link{RunNMF}} function. The structure of WH.list should be:
#'
#' \itemize{
#'  \item\code{sample_1}
#'  \itemize{
#'  \item\code{W} Normalized W matrix from NMF
#'  \item\code{H} Normalized H matrix from NMF
#' }
#' \item\code{sample_2}
#'  \itemize{
#'  \item\code{W} Normalized W matrix from NMF
#'  \item\code{H} Normalized H matrix from NMF
#' }
#' ...
#' \item\code{sample_n}
#' }
#'
#' @param cluster.result The result of \code{\link{OvlpHeatmap}} function contains cluster reults for each program. Or a named vector of cluster result.
#' @param top The number of top genes to extract for each meta-program. Default is 50.
#' @param key A prefix string for the names of extracted meta-programs. Default is 'MP_'.
#' @param only.gene Logical, if TRUE, only gene names are returned; if FALSE, named gene scores are returned. Default is TRUE.
#'
#' @return A list of extracted meta-programs.
#' Each element of the list represents a meta-program and contains either gene names or gene scores (if only.gene = FALSE).
#'
#'
#' @importFrom utils head
#' @details
#' The function extracts meta-programs from the NMF results by aggregating the results from multiple samples.
#' It first aligns the genes across all samples and then computes the average expression profile for each meta-program.
#' The top genes with the highest average expression scores are selected and returned.
#'
#' @seealso
#' \code{\link{RunNMF}} for the function to perform NMF analysis.
#'
#' @export
#'

MetaProgram = function(WH.list, cluster.result, top = 50, key = 'MP_', only.gene = TRUE ){
    #extra the cluster result
    if(is.list(cluster.result)){
        cluster = cluster.result$Cluster
    }else{
        cluster = cluster.result
    }

    WH.list = WH.list[!sapply(WH.list,is.null)]
    all_gene = lapply(WH.list, function(WH){
        #those are hvg in an individual
        rownames(WH$W)
    }) %>% Reduce(f = union)

    all_pg = lapply(WH.list, function(WH){
        colnames(WH$W)
    }) %>% Reduce(f = union)

    if(!all(names(cluster) %in% all_pg)){
        warning('Can not find program: ', paste(setdiff(names(cluster), all_pg), collapse = ' '), ' in the WH.list')
    }
    all_pg = intersect(all_pg, names(cluster))
    cluster = cluster[all_pg]

    all_W = lapply(WH.list, function(WH){
        neo_W = matrix(0, ncol = ncol(WH$W), nrow = length(all_gene))
        rownames(neo_W) = all_gene
        colnames(neo_W) = colnames(WH$W)
        neo_W[match(rownames(WH$W), all_gene),] = WH$W
        neo_W = neo_W[,intersect(all_pg, colnames(neo_W)), drop = FALSE]
    }) %>% do.call(what = cbind)

    ls_mp = split(names(cluster), cluster) %>% lapply(function(pgs){
        mat_pgs = all_W[,pgs]
        genes = rowMeans(mat_pgs) %>% sort(decreasing = TRUE) %>% head(n = top)
        if(only.gene){
            genes = names(genes)
        }

        return(genes)
    })

    names(ls_mp) = paste0(key, names(ls_mp))
    return(ls_mp)
}
