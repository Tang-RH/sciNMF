#' Cluster Metrics Plot
#'
#' This function is used to plot the clustering results of provided different numbers of clusters,
#' It calculates and visualizes key clustering evaluation metrics to aid in the assessment of clustering performance.
#'
#' @param mat.ovlp A symmetric numeric matrix generated by \code{\link{OverlapMat}}.
#' @param num.clusters The number range of clusters used for clustering analysis. The function will
#' calculate evaluation metrics for each specified number of clusters. Default is 2:20.
#' @param ncol The number of columns for arranging subplots in the plot area. Default is 3.
#' @param distance.clustering The distance measure used for clustering. Possible values are \code{"Intersection"} for intersection between two programs as similarity,
#' \code{"Jaccard"} for Jaccard similarity, \code{"correlation"} for Pearson correlation and all the distances supported by \code{\link[stats]{dist}}, such as \code{"euclidean"}, etc.
#' Default is \code{"Intersection"}. If this parameter is set as \code{"Jaccard"}, the \code{Jaccard} parameter in \code{\link{OverlapMat}} should be set as \code{TRUE}.
#' @param max.intersect When the \code{"distance.clustering"} is set to \code{"Intersection"}, the distance between two programs will be calculated as the difference of \code{max.intersect} and \code{intersection}. Default is 50.
#' @param method.clustering The clustering method. Possible values are all the methods supported by \code{\link[stats]{hclust}},
#'  such as \code{"ward.D2"}, \code{"single"}, \code{"average"}, etc. Default is \code{"ward.D2"}.
#'
#' @return Returns a plot displaying cluster metrics.
#'  The plot includes lines and scatter points representing clustering evaluation metrics for different numbers of clusters.
#'
#' @importFrom ggplot2 ggplot geom_line geom_point theme_bw facet_wrap
#' @importFrom fpc cqcluster.stats
#' @importFrom stats hclust dist as.dist cutree
#'
#' @examples
#'
#' # Sample transcriptional programs
#' set.seed(123)
#' ls_pg <- lapply(1:50, function(i) {
#'   paste0("g", sample(1:1000, 50))
#' })
#'
#' # Generate the ovelap matrix
#' ovlp <- OverlapMat(ls_pg)
#' # Generate cluster plot
#' cluster_plot <- ClusterMetricsPlot(ovlp, num.clusters = 2:10, ncol = 2)
#'
#' # Display the cluster plot
#' print(cluster_plot)
#'
#' @seealso
#' \code{\link{OverlapMat}}, \code{\link[stats]{hclust}}, \code{\link[stats]{dist}},  \code{\link[fpc]{cqcluster.stats}}
#'
#' @export
#'
ClusterMetricsPlot <- function(mat.ovlp, num.clusters = 2:20, ncol = 3,
                               distance.clustering = "Intersection", max.intersect = 50,
                               method.clustering = "ward.D2") {
  if (distance.clustering == "Intersection") {
    res_dist <- as.dist(max.intersect - mat.ovlp)
  } else if (distance.clustering == "Jaccard") {
    res_dist <- as.dist(1 - mat.ovlp)
  } else if (distance.clustering == "correlation") {
    res_dist <- as.dist(1 - cor(mat.ovlp))
  } else {
    res_dist <- dist(mat.ovlp, method = distance.clustering)
  }

  res_cluster <- hclust(res_dist, method = method.clustering)
  # avoid invalid cluster number
  max_cluster_number <- nrow(mat.ovlp) - 1

  num.clusters <- intersect(num.clusters, 2:max_cluster_number)

  df_pl <- lapply(num.clusters, function(k) {
    cluster <- stats::cutree(res_cluster, k)
    res_evaluation <- fpc::cqcluster.stats(res_dist, cluster)
    sub_df <- data.frame(
      Index = c(
        "SilhouetteWidth", "Calinski-Harabasz", "SeparationIndex", "WidestGap",
        "PearsonGamma", "DunnIndex"
      ),
      Value = c(
        res_evaluation$asw, res_evaluation$ch, res_evaluation$sindex, res_evaluation$widestgap,
        res_evaluation$pearsongamma, res_evaluation$dunn
      ),
      ClusterNumber = k
    )
  }) %>% do.call(what = rbind)

  df_pl$Index <- factor(df_pl$Index, levels = unique(df_pl$Index))
  pl <- ggplot(df_pl, aes(ClusterNumber, Value)) +
    geom_line(col = "blue") +
    geom_point(col = "red") +
    theme_bw() +
    facet_wrap(~Index, ncol = ncol, scales = "free")
  return(pl)
}
