#' Title
#'
#' @param x
#' @param k
#' @param seq.size
#'
#' @return
#'
#' @examples
ExtractSubSequence <- function(x=NULL,k = 20,seq.size = 201){
    extend <- ((seq.size-1)/2)
    center <- seq(1,length(x),k)
    start <- center - extend
    end <- center + extend
    Viewseq <- Biostrings::Views(x, start=start[start>0&end<=length(x)], end=end[start>0&end<=length(x)])
    sequences <- Biostrings::DNAStringSet(Viewseq)

    results <- cbind(as.data.frame(IRanges::ranges(Viewseq)),seq=as.character(Viewseq))
    return(results)
}
