#' Cluster Programs Based on Overlap Matrix
#'
#' This function clusters robust programs using the provided overlap matrix
#' and filters the clusters based on their size (number of robust programs).
#'
#' @param mat.ovlp A matrix representing the overlap between programs.
#' @param cut.num Number of clusters to be created.
#' @param min.cluster.size Minimum size of a cluster to be retained default is 0.
#' @param method.clustering Method for hierarchical clustering.
#'   Options include "ward.D2" (default), "single", "complete" and other methods supported by \code{\link[stats]{hclust}}
#' @param distance.clustering Method for calculating distances between clusters.
#'   Options include "Intersection" (default), "Jaccard", "correlation", and any other valid method for \code{\link[stats]{dist}}.
#' @param max.intersect Maximum intersection value (used when distance.clustering = "Intersection"). Default is 50.
#'
#' @return A vector indicating the cluster assignment for each program.
#'
#' @examples
#' # Example lists
#' set.seed(123)
#' ls_pg <- lapply(setNames(nm = paste0("PG", 1:50)), function(i) {
#'   paste0("g", sample(1:100, 20))
#' })
#'
#' # Calculate the overlap matrix
#' ovlp <- OverlapMat(ls_pg)
#' clusters <- ClusterPG(ovlp, cut.num = 5, min.cluster.size = 4, max.intersect = 20)
#' print(clusters)
#'
#' @seealso
#' \code{\link[stats]{hclust}}, \code{\link[stats]{cutree}}, \code{\link[stats]{dist}}
#'
#' @import stats
#'
#' @export
#'
ClusterPG <- function(mat.ovlp, cut.num, min.cluster.size = 0,
                      method.clustering = "ward.D2",
                      distance.clustering = "Intersection", max.intersect = 50) {
  if (distance.clustering == "Intersection") {
    res_dist <- as.dist(max.intersect - mat.ovlp)
  } else if (distance.clustering == "Jaccard") {
    res_dist <- as.dist(1 - mat.ovlp)
  } else if (distance.clustering == "correlation") {
    res_dist <- as.dist(1 - cor(mat.ovlp))
  } else {
    res_dist <- dist(mat.ovlp, method = distance.clustering)
  }

  res_hclust <- hclust(res_dist, method = method.clustering)
  cluster <- cutree(res_hclust, cut.num)[res_hclust$order]
  df_freq <- data.frame(table(cluster, dnn = "Cluster"))
  df_freq <- subset(df_freq, Freq >= min.cluster.size)
  cluster <- cluster[cluster %in% df_freq$Cluster]
  return(cluster)
}
