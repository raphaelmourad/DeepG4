#' ExtractSubSequencefromBed : Internal function for DeepG4Scan.
#'
#' @param x An object of class GRanges.
#' @param x.atac A SimpleRleList object containing DNA accessibility.
#' @param k size of the sliding windows.
#' @param GENOME a BSgenome object containing  the DNA sequence of genome of interest.
#' @param seq.size numeric value representing the sequence size accepted by our model. Don't change it unless you want to use our function with a custom model.
#' @param nb.threads number of threads to use, default 1.
#' @param use.bg a boolean. Set to \code{TRUE} you want to normalize the accessibility using a windows background of windows_bg.
#' @param windows_bg numeric value who define the windows use to get background signal.
#' @param treshold_bg numeric value who set the treshold signal/background.
#'
#' @return A data.frame with position of subsequence, DNA sequence and accessibility.
#' @export
#'
#' @examples
ExtractSubSequencefromBed <- function(x=NULL,x.atac=NULL,k = 20,GENOME=NULL,seq.size = 201,nb.threads=1,use.bg=T,windows_bg=5000,treshold_bg = 2){
    extend <- ((seq.size-1)/2)
    windows = IRanges::slidingWindows(x, width=k,step = k)
    res <- parallel::mclapply(1:length(windows),function(i){

        swindows <- windows[[i]]
        end(swindows) <- BiocGenerics::start(swindows) + extend
        start(swindows) <- BiocGenerics::start(swindows) - extend

        cov <- x.atac[[as.character(GenomeInfoDb::seqnames(x[i]))]]

        score <- IRanges::Views(cov, start = BiocGenerics::start(swindows), end = BiocGenerics::end(swindows))
        score <- rowMeans(as.matrix(score))
        score[is.na(score)] <- 0

        if(use.bg){
            swindows.bg <- resize(swindows,windows_bg,fix="center")
            score.bg <- IRanges::Views(cov, start = BiocGenerics::start(swindows.bg), end = BiocGenerics::end(swindows.bg))
            score.bg <- rowMeans(as.matrix(score.bg))
            score.bg[is.na(score.bg)] <- 0
            my_test <- (score/score.bg)<treshold_bg
            score[my_test] <- 0
        }
        results <- cbind(seqnames = as.character(seqnames(x[i])),as.data.frame(IRanges::ranges(swindows)),score=score,seq=as.character(Biostrings::getSeq(GENOME,swindows)))

        return(results)
    },mc.cores=nb.threads)

    return(do.call("rbind",res))
}
