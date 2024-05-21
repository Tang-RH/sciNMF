#' Robust Programs Identification
#'
#' This function identifies robust gene programs using a modified method proposed in the Nature 2023 paper.
#' The method aims to remove redundant programs that are highly similar within the same patient.
#'
#' @param WH.list
#' A nested list containing NMF results for each individual, generated from the \code{\link{RunNMF}} function. The structure of WH.list should be:
#'
#' \itemize{
#'  \item\code{sample_1}
#'  \itemize{
#'  \item\code{W} W matrix from NMF
#'  \item\code{H} H matrix from NMF
#' }
#' \item\code{sample_2}
#'  \itemize{
#'  \item\code{W} W matrix from NMF
#'  \item\code{H} H matrix from NMF
#' }
#' ...
#' \item\code{sample_n}
#' }
#'
#' @param top The number of top genes to extract for each robust program. Default is 50.
#' @param IQR.cut The threshold for the interquartile range (IQR) of gene usage to filter programs during QC. Default is 0.1.
#' @param median.cut The threshold for the median gene usage to filter programs during QC. Default is 0.02.
#' @param intra.min The threshold to identify similar programs within an individual
#' Number of overlap genes greater than or equal to \code{intra.min} are considered similar within the same individual
#' Default is 35
#' @param intra.rep The minimum number of similar programs. If a program has less than \code{intra.rep} similar programs in an individual,
#' it will be removed. Default is 1
#' @param inter.filter Logical value indicates whether filter the programs based on inter-tumor similarities. Default is TRUE
#' @param inter.min The threshold to identify robust programs across individuals
#' Number of overlap genes greater than or equal to \code{inter.min} are considered similar across individuals
#' Default is 10
#' @param inter.rep The minimum number of similar programs across other individuals to be considered as a robust program after intra-tumor filtering.
#' Default is 1
#' @param intra.max The threshold to identify redundant programs within an individual.
#' Number of overlap genes greater than or equal to \code{intra.max} are considered redundant within the same individual.
#' Default is 10
#'
#' @return A list containing robust gene programs for each individual
#' Each element of the list represents a program and contains \code{top} genes with gene names.
#'
#' @details
#' The function performs Quality Control (QC) on the programs based on the Interquartile Range (IQR) and median usage (the intra-rank-normalized H matrix).
#' Programs with low IQR and median usage are filtered out.
#'
#' Following QC, it removes programs that have fewer than \code{intra.rep} similar programs (number of overlapping genes between two programs >= \code{intra.min}) within the same individual (for selecting robust programs).
#'
#' Subsequently, it filters out programs which have fewer than \code{inter.rep} similar programs (number of overlapping genes between two programs >= \code{inter.min}) across individuals.
#'
#' Finally, it further filters out redundant programs, selecting only programs that have an intersection smaller than \code{intra.max} with previously selected programs within an individual.
#'
#' The remaining programs are considered robust and are returned as the final result.
#'
#'
#' @seealso
#' \code{\link{RunNMF}} for the function to perform NMF analysis.
#'
#' @importFrom stats setNames quantile
#'
#' @references
#' Gavish, Avishai et al. “Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours.” Nature vol. 618,7965 (2023): 598-606. \url{https://doi.org/10.1038/s41586-023-06130-4}
#' @export
#'


RobustProgram = function(WH.list, top = 50, IQR.cut = 0.1, median.cut = 0, intra.min = 35, intra.rep = 1,
                         inter.filter = TRUE,inter.min = 10, inter.rep = 1, intra.max = 10){
    
    WH.list = WH.list[!sapply(WH.list, is.null)]
    #QC by IQR and median useage
    ls_pat_pg = lapply(WH.list, function(WH){
        #normalize H by col to calculate IQR for each k
        Ks = sub(".*_K([0-9]+)_P[0-9]+$", "\\1", rownames(WH$H))
        H = split(data.frame(WH$H, check.names = FALSE), Ks) %>% lapply(function(sub_H){
            H_ratio = apply(sub_H, 2, function(me){
                me/sum(me)
            })
            return(H_ratio)
        }) %>% {names(.) = NULL
        do.call(what = rbind, args = .)}
        
        mat_quat = apply(H,1, quantile)
        idx_median = mat_quat['50%',] >= median.cut
        idx_IQR = mat_quat['75%',] - mat_quat['25%',] >= IQR.cut
        
        W_filter = WH$W[, idx_median & idx_IQR]
        
        pgs = lapply(setNames(nm = colnames(W_filter)), function(pg){
            pg = W_filter[,pg]
            head(sort(pg, decreasing = TRUE), n = top)
        })
        #filter by intra.min
        mat_ovlp = OverlapMat(lapply(pgs, names))
        #the mat_ovlp contains the program itself, so > intra.rep
        idx_pg = apply(mat_ovlp, 1, function(pgop){sum(pgop >= intra.min)}) > intra.rep
        if(sum(idx_pg) == 0){
            #some samples can not pass the intra.min >= intra.rep
            return(NULL)
        }
        pgs = pgs[idx_pg]
    })
    
    ls_pat_pg = ls_pat_pg[!sapply(ls_pat_pg, is.null)]
    
    #filter inter.min
    mat_ovlp_all = OverlapMat(lapply(unlist(ls_pat_pg, recursive = FALSE), names))
    ls_pat_names = lapply(ls_pat_pg, names)
    colnames(mat_ovlp_all) = rownames(mat_ovlp_all) = unlist(ls_pat_names, use.names = FALSE)
    
    ls_keep_pg = lapply(1:length(ls_pat_pg), function(pat){
        idx_pgs = match(names(ls_pat_pg[[pat]]), colnames(mat_ovlp_all))
        sub_mat_ovlp = mat_ovlp_all[idx_pgs, -idx_pgs, drop = FALSE]
        #the sub_mat_ovlp doesn't contain the program itself, so >= inter.rep
        idx_inter = apply(sub_mat_ovlp, 1, function(pgop){ sum(pgop >= inter.min) >= inter.rep })
        
        #filter out without inter tumor rep
        pgs = names(idx_inter)[idx_inter]
        #some samples can not pass the inter.min >= inter.rep
        if(length(pgs) == 0){
            return(NULL)
        }
        pgs = sort(apply(sub_mat_ovlp[pgs, , drop = FALSE], 1, max), decreasing = TRUE)
        
        if(length(pgs) > 1){
            keep_pgs = names(pgs)[1]
            for(pg_test in names(pgs)[-1]){
                if(max(mat_ovlp_all[keep_pgs,pg_test]) <= intra.max ){
                    keep_pgs = c(keep_pgs, pg_test)
                }
            }
        }else{
            keep_pgs = names(pgs)
        }
        
        return(lapply(ls_pat_pg[[pat]][keep_pgs], names))
        
    }) %>% unlist(recursive = FALSE)
    
    return(ls_keep_pg[!sapply(ls_keep_pg, is.null)])
}
