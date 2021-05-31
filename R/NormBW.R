#' Title
#'
#' @param x
#' @param binbed
#'
#' @return
#' @export
#'
#' @examples
NormBW <- function(x,binbed){
    ranges_bin <-  base::range(IRanges::subsetByOverlaps(x,binbed)$score,na.rm = TRUE, finite = TRUE)
    x$score <- scales::rescale(x$score,to = c(0,1),from = ranges_bin)
    x <- IRanges::coverage(x,weight = "score")
    return(x)
}
