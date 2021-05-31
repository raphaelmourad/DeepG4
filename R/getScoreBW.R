#' Title
#'
#' @param WIG
#' @param BED
#'
#' @return
#' @export
#'
#' @examples
getScoreBW <- function (WIG, BED)
{
    res <- lapply(split(BED, droplevels(GenomeInfoDb::seqnames(BED))), function(zz) {
        cov <- WIG[[unique(as.character(GenomeInfoDb::seqnames(zz)))]]
        score <- IRanges::Views(cov, start = BiocGenerics::start(zz), end = BiocGenerics::end(zz))
        return(IRanges::viewMeans(score))
    })
    res <- do.call("c",res)
    return(res)
}
