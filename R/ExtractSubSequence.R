#' ExtractSubSequencefromBed : Internal function for DeepG4Scan.
#'
#' @param x An object of class GRanges.
#' @param k size of the sliding windows.
#' @param seq.size numeric value representing the sequence size accepted by our model. Don't change it unless you want to use our function with a custom model.
#'
#' @return A data.frame with position of subsequence and DNA sequence.
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
