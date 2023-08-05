#' IQRPlot: Generate Boxplots with IQR and Median Filter
#'
#' The IQRPlot function is used to generate boxplots for each program based on the IQR (Interquartile Range) and median filter. It calculates the IQR and median usage of each program from the NMF (Non-Negative Matrix Factorization) results and filters out programs with low IQR and median values. The higher IQR indicates a larger range of usage, which is useful for quality control.
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
#' @param IQR.cut The cutoff value for IQR. Programs with IQR above this value will be kept in the plot.
#' @param median.cut The cutoff value for median. Programs with median above this value will be kept in the plot.
#' @param grid A logical value indicating whether to arrange the plots of samples in a grid (TRUE) or return a list of ggplot objects (FALSE).
#' @param ncol Determines the number of columns when arranging the plots in a grid. If NULL, it will be automatically determined based on the number of samples.
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
IQRPlot = function(WH.list, IQR.cut = 0.1, median.cut = 0, grid = TRUE, ncol = NULL){
    WH.list = WH.list[!sapply(WH.list, is.null)]
    ls_pl = lapply(WH.list, function(WH){

        df_pl = tidyr::gather(data.frame(t(WH$H)),'Program','Ratio')
        df_pl$Program = rep(rownames(WH$H), each = ncol(WH$H))
        mat_quat = apply(WH$H,1, quantile)
        idx_median = mat_quat['50%',] > median.cut
        idx_IQR = mat_quat['75%',] - mat_quat['25%',] > IQR.cut
        df_pl$Keep = ifelse(df_pl$Program %in% colnames(mat_quat)[idx_median & idx_IQR], 'TRUE', 'FALSE')
        df_pl$K = sub(".*_K([0-9]+)_P[0-9]+$", "\\1", df_pl$Program)
        sam = gsub("_K[0-9]+_P[0-9]+$", "", df_pl$Program[1])
        df_pl$Program = factor(df_pl$Program, levels = colnames(mat_quat)[order(mat_quat['50%',])])

        p = ggplot(df_pl, aes(Program, Ratio)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(aes(color = Keep), alpha=0.4,width = 0.2, size = 0) +
          theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) +
          ggtitle(paste0(sam, '_k',min(as.numeric(df_pl$K)) ,'to', max(as.numeric(df_pl$K)))) +
          facet_grid(~K,drop = TRUE,scales = "free",space = 'free') +
          theme(plot.title = element_text(hjust = 0.5))

        return(p)

    })
    if(grid){
        if(is.null(ncol)){ncol = round(sqrt(length(WH.list)))}
        ls_pl = cowplot::plot_grid(plotlist = ls_pl, ncol = ncol)
    }
    return(ls_pl)
}
