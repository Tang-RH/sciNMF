#' IQRPlot: Generate Boxplots with IQR and Median Filter
#'
#' The IQRPlot function is used to generate boxplots for each program based on the IQR (Interquartile Range) and median filter. 
#' The IQR and median usage of each program is calculated from normalized H matrices, the NMF (Non-Negative Matrix Factorization) results. 
#' The higher IQR indicates a larger range of usage, which is useful for quality control.
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
#' @param IQR.cut The cutoff value for IQR. Programs with IQR above this value will be kept in the plot.
#' @param median.cut The cutoff value for median. Programs with median above this value will be kept in the plot.
#' @param grid A logical value indicating whether to arrange the plots of samples in a grid (TRUE) or return a list of ggplot objects (FALSE).
#' @param ncol Determines the number of columns when arranging the plots in a grid. If NULL, it will be automatically determined based on the number of samples.
#' @param align (Optional) Specifies whether graphs in the grid should be horizontally ("h") or vertically ("v") aligned. Options are "none" (default), "hv" (align in both directions), "h", and "v".
#' @details
#' First, all the program usages for each cell will be normalized to 1. Then the IQR and median usage of each program will be calculated for quality control.
#' The programs with IQR and median usage below the cutoff will be labeled as 'Remove'
#' 
#' @return A girded ggplot object or a list of ggplot objects containing boxplot for each program with IQR and median filtering.
#'
#' @import ggplot2
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom cowplot plot_grid
#' @importFrom stats quantile
#' @importFrom utils tail
#' @examples
#' \dontrun{
#'     result <- IQRPlot(WH.list, IQR.cut = 0.1, median.cut = 0, grid = TRUE)
#' }
#'
#' @keywords NMF result boxplot visualization
#'
#' @export
#'
IQRPlot = function(WH.list, IQR.cut = 0.1, median.cut = 0.02, grid = TRUE, ncol = NULL, align = c("hv",  "h", "v", "none")){
    WH.list = WH.list[!sapply(WH.list, is.null)]
    ls_pl = lapply(WH.list, function(WH){
        #normalize H by col to calculate IQR for each k
        Ks = sub(".*_K([0-9]+)_P[0-9]+$", "\\1", rownames(WH$H))
        H = split(data.frame(WH$H, check.names = FALSE), Ks) %>% lapply(function(sub_H){
            H_ratio = apply(sub_H, 2, function(me){
                me/sum(me)
            })
        }) %>% {names(.) = NULL
        do.call(what = rbind, args = .)}
        

        df_pl = tidyr::gather(data.frame(t(H)),'Program','Ratio')
        mat_quat = apply(H,1, quantile)
        idx_median = mat_quat['50%',] > median.cut
        idx_IQR = mat_quat['75%',] - mat_quat['25%',] > IQR.cut
        df_pl$Keep = ifelse(df_pl$Program %in% colnames(mat_quat)[idx_median & idx_IQR], 'Keep', 'Remove')
        df_pl$K = sub(".*_K([0-9]+)_P[0-9]+$", "\\1", df_pl$Program)
        sam = gsub("_K[0-9]+_P[0-9]+$", "", df_pl$Program[1])
        df_pl$Program = factor(df_pl$Program, levels = colnames(mat_quat)[order(mat_quat['50%',])])
        col_box = c('Keep' = '#1F77B4', 'Remove' = '#D62728')

        p = ggplot(df_pl, aes(Program, Ratio)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(aes(color = Keep), alpha=0.4,width = 0.2, size = 0) +
          theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) +
          ggtitle(paste0(sam, '_k',min(as.numeric(df_pl$K)) ,'to', max(as.numeric(df_pl$K)))) +
          facet_grid(~K,drop = TRUE,scales = "free",space = 'free') +
          theme(plot.title = element_text(hjust = 0.5)) + 
          scale_color_manual(values = col_box)

        return(p)

    })
    if(grid){
        if(is.null(ncol)){ncol = round(sqrt(length(WH.list)))}
        ls_pl = cowplot::plot_grid(plotlist = ls_pl, ncol = ncol, align = align)
    }
    return(ls_pl)
}
