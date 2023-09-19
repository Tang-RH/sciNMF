#' Calculate Overlap Matrix
#'
#' The \code{OverlapMat} function calculates the overlap between two lists of gene sets and creates an overlap matrix.
#'
#' @param ls1 The first list of gene sets for which the overlap is to be calculated.
#' @param ls2 The second list of gene sets for which the overlap is to be calculated. By default, it is set to \code{ls1}, which means the function calculates the overlap of \code{ls1} with itself.
#' @param Dice A logical value indicating whether to calculate Dice similarity. Default is FALSE.
#' @param Jaccard A logical value indicating whether to calculate Jaccard similarity. Default is FALSE.
#'
#' @details
#' The \code{OverlapMat} function takes two lists of gene sets and creates a matrix with rows representing elements from \code{ls1} and columns representing gene sets from \code{ls2}. The function calculates the size of the intersection between each pair of gene sets from \code{ls1} and \code{ls2}. If \code{Jaccard} is \code{TRUE}, it also calculates the Jaccard similarity by dividing the size of the intersection by the size of the union of the two sets.
#'
#' @return A matrix representing the overlap between gene sets in the input lists.
#'
#' @examples
#' # Example lists
#' set.seed(123)
#' ls_pg <- lapply(1:50, function(i){
#'   paste0('g',sample(1:100,20))
#' })
#'
#' # Calculate the overlap matrix
#' ovlp <- OverlapMat(ls_pg)
#'
#' # Print the result
#' print(ovlp)
#'
#' @export
#'
OverlapMat = function(ls1, ls2 = ls1, Dice = FALSE, Jaccard = FALSE ){
    if(Dice & Jaccard){warning('Both Dice and Jaccard are TRUE, return Dice')}
    mat = sapply(ls2, function(m2) {
        sapply(ls1, function(m1) {
            inset = intersect(m1, m2)
            if(Dice){
                res = length(inset)*2/length(c(m1,m2))
            }else if(Jaccard){
                uniset = union(m1, m2)
                res = length(inset)/length(uniset)
            }else{
                res = length(inset)
            }
            return(res)
        })
    })
    return(mat)
}
