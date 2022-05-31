#' Pre-Computing function of DeepG4. To use before DeepG4 main function
#'
#' @param BED An object of class GRanges.
#' @param ATAC A character path of bigWig/bedGraph file, or an object of class GRanges/SimpleRleList.
#' @param is.bw a boolean. Set to \code{TRUE} if you want to use rtracklayer::import.bw fonction, \code{FALSE} use rtracklayer::import.bedGraph.
#' @param GENOME a BSgenome object containing  the DNA sequence of genome of interest.
#' @param use.bg a boolean. Set to \code{TRUE} you want to normalize the accessibility using a windows background of windows_bg.
#' @param windows_bg numeric value who define the windows use to get background signal.
#' @param treshold_bg numeric value who set the treshold signal/background.
#'
#' @return a list(DNAStringSet,vector) with fasta sequence from GENOME of BED, and vector of accessibility of BED.
#' @export
#'
#' @examples
DeepG4InputFromBED <- function(BED = NULL,ATAC = NULL,is.bw = TRUE,GENOME = NULL,use.bg = TRUE,windows_bg=5000,treshold_bg = 2){
    if (is.null(GENOME)) {
        stop("GENOME must be provided (see ?DeePG4InputFromBED for accepted formats).",
             call. = FALSE)
    }
    if (class(GENOME) != "BSgenome") {
        stop("GENOME must be a BSgenome object.",
             call. = FALSE)
    }
    if (is.null(BED)) {
        stop("BED must be provided (see ?DeePG4InputFromBED for accepted formats).",
             call. = FALSE)
    }
    if(class(BED) == "character"){
        BED <- tryCatch(
            {
                rtracklayer::import.bed(BED)
            },
            error=function(cond) {
                message(paste(BED," is not in bed format."))
                message(cond)
                # Choose a return value in case of error
                return(NULL)
            })
        if (is.null(BED)) {
            stop("BED must be provided in a valid format (see ?DeePG4InputFromBED for accepted formats).",
                 call. = FALSE)
        }
    }else if(class(BED) != "GRanges"){
        stop("BED must be provided in a valid format (see ?DeePG4InputFromBED for accepted formats).",
             call. = FALSE)
    }
    #ATAC TESTING
    if (is.null(ATAC)) {
        stop("ATAC must be provided (see ?DeePG4InputFromBED for accepted formats).",
             call. = FALSE)
    }
    if(class(ATAC) == "character"){
        if(is.bw){
            ATAC <- tryCatch(
                {
                    rtracklayer::import.bw(ATAC)
                },
                error=function(cond) {
                    message(paste(ATAC," is not in BigWig format.\nIf you want to load a bedGraph file, set is.bw to FALSE:"))
                    message(cond)
                    # Choose a return value in case of error
                    return(NULL)
                })
            if (is.null(ATAC)) {
                stop("ATAC must be provided in a valid format (see ?DeePG4InputFromBED for accepted formats).",
                     call. = FALSE)
            }
        }else{
            ATAC <- tryCatch(
                {
                    rtracklayer::import.bedGraph(ATAC)
                },
                error=function(cond) {
                    message(paste(ATAC," is not in bedGraph format.\nIf you want to load a bigwig file, set is.bw to TRUE:"))
                    message(cond)
                    # Choose a return value in case of error
                    return(NULL)
                })
            if (is.null(ATAC)) {
                stop("ATAC must be provided in a valid format (see ?DeePG4InputFromBED for accepted formats).",
                     call. = FALSE)
            }
        }
    }else if(class(ATAC) == "SimpleRleList"){
        ATAC <- as(ATAC,"GRanges")
    }else if(class(ATAC) == "GRanges"){
        if(!"score" %in% colnames(values(ATAC))){
            stop("ATAC must be a GRanges object with a metadata column 'score'",
                 call. = FALSE)
        }
        if(class(ATAC$score) != "numeric"){
            stop("score column must be a numeric vector.",
                 call. = FALSE)
        }
    }else{
        stop("ATAC must be a character (path for a file in BedGraph/BigWig format), or a SimpleRleList/GRanges class",
             call. = FALSE)
    }
    #Normalize ATAC-seq using the previously computed bins
    BED$order <- 1:length(BED)
    X <- Biostrings::getSeq(GENOME,BED)
    BED <- BiocGenerics::sort(GenomeInfoDb::sortSeqlevels(BED))

    binbed <- rtracklayer::import.bed(system.file("extdata", "random_region_for_scaling_min_max.bed", package = "DeepG4"))
    ATAC <- NormBW(ATAC,binbed)
    X.ATAC <- getScoreBW(ATAC,BED)
    X.ATAC[is.na(X.ATAC)] <- 0

    if(use.bg){
        BED.bg <- resize(BED,windows_bg,fix="center")
        X.ATAC.bg <- getScoreBW(ATAC,BED.bg)
        X.ATAC.bg[is.na(X.ATAC.bg)] <- 0
        my_test <- (X.ATAC/X.ATAC.bg)<treshold_bg
        X.ATAC[my_test] <- 0
    }

    X.ATAC <- X.ATAC[order(BED$order)]

    return(list(X,as.vector(X.ATAC)))
}
