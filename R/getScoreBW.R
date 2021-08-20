#' getScoreBW extract DNA accessibility from WIG to BED positions.
#'
#' @param WIG An object of class GRanges containing DNA accessibility.
#' @param BED An object of class GRanges.
#'
#' @return A SimpleRleList object containing DNA accessibility.
#' @export A numeric vector of DNA accessibility
#'
#' @examples
getScoreBW <- function (WIG, BED,forScan=F)
{
    res <- do.call("rbind",lapply(split(BED, droplevels(GenomeInfoDb::seqnames(BED))), function(zz) {
        cov <- WIG[[unique(as.character(GenomeInfoDb::seqnames(zz)))]]
        score <- IRanges::Views(cov, start = BiocGenerics::start(zz), end = BiocGenerics::end(zz))
        return(as.matrix(score))
    }))
    if(forScan){
        return(res)
    }else{
        return(rowMeans(res))
    }
}
