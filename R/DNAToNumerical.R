#' Convert DNAStringSet object into one-hot encoding
#'
#' @param x DNAStringSet object
#' @param tabv named vector list of numerical values which indicate the numerical value of a nucleotide.
#' @param lower.case Set to \code{TRUE} if the DNA in your sequence is encoded in lower.case i.e. "acgt", default set to \code{FALSE}.
#' @param seq.size Set the sequence maximal length value authorized by our model (default to 201).
#'
#' @return An array of dimension \code{nrow(x),ncol(x),length(tabv)}
#' @export
#' @examples
#' x <- Biostrings::DNAStringSet(c("ACGT"))
#' x_onehot <- DNAToNumerical(x)
#' x_onehot
DNAToNumerical <- function(x,tabv = c("N"=5,"T"=4,"G"=3,"C"=2,"A"=1),lower.case=F,seq.size = 201){
    if(lower.case){
        names(tabv) <- tolower(tabv)
    }
    x <- Biostrings::as.matrix(x)
    listMat <- list()
    for(i in 1:length(tabv)){
        nuc_index <- tabv[[i]]
        nuc_value <- names(tabv[i])
        mat <- matrix(0,nrow(x),ncol(x))
        mat[x==nuc_value] <- 1
        if(ncol(x)<seq.size){
            mat <- cbind(mat,matrix(0,nrow(x),seq.size-ncol(x)))
        }
        listMat[[nuc_index]] <- mat
    }
    arrayout <- array(unlist(listMat), dim = c(nrow(listMat[[1]]), ncol(listMat[[1]]), length(listMat)))
    return(arrayout)

}
