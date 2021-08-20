#' Normalize a coverage file (from bigWig or bedGraph) using scales::rescale from a specific set of ranges to set values between [0,1]
#'
#' @param x An object of class GRanges, the coverage file to be normalized
#' @param binbed An object of class GRanges, the genomic position used to normalize x.
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
