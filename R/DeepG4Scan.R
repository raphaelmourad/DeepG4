#' Scanning of potential active G4 on a sequence or a list of sequences using a sliding window of size k.
#'
#' @param X An object of class character,list or DNAStringSet/DNAStringSetList with DNA sequences.
#' @param k size of the sliding windows.
#' @param treshold  numeric value who define the treshold to use to consider a sequence asc ontaining an active G4.
#' @param threads  numeric value who define the number of threads used in DeepG4Scan (Generate sub sequences)
#' @param lower.case a boolean. Set to \code{TRUE} if elements of X are in lower case (default to FALSE).
#'
#' @return a data.frame with the position of potential active G4 across input sequences.
.callDeepG4Scan <- function(X=NULL,k=20,treshold = 0.5,threads = 1,lower.case=F){
    seq.size = 201
    # Packages check
    if (!requireNamespace("keras", quietly = TRUE)) {
        stop("Package \"keras\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (!requireNamespace("Biostrings", quietly = TRUE)) {
        stop("Package \"Biostrings\" needed for this function to work. Please install it.",
             call. = FALSE)
    }

    ## Check model class and convert into DNAStringSet object if needed
    if(!class(X)[[1]] %in%c("DNAString","DNAStringSet","DNAStringSetList")){
        if(class(X) == "character"){
            X <- Biostrings::DNAStringSet(X)
        }else if(class(X) == "list"){
            if(class(X[[1]])[[1]] == "DNAString"){
                X <- as(X,"DNAStringSet")
            }else if(class(X[[1]])[[1]] == "character"){
                X <- Biostrings::DNAStringSet(unlist(X))
            }else{
                stop("X must be a list of DNAString/character class",
                     call. = FALSE)
            }
        }else{
            stop("X must be a character, a list or a DNAStringSet/DNAStringSetList class",
                 call. = FALSE)
        }
    }else if(class(X)[[1]] =="DNAStringSetList"){
        X <- unlist(Biostrings::DNAStringSetList(X))
    }
    ## Check DNA composition
    if(class(X)[[1]] !="DNAString"&&length(X)>1){
        results <- parallel::mclapply(1:length(X),function(i){
            x <- X[i][[1]]
            results <- cbind(seqnames= i,ExtractSubSequence(x=x,k=k,seq.size = seq.size))
            return(results)
        },mc.cores = threads)
        results <- do.call(rbind,results)
    }else{
        if(class(X)[[1]] =="DNAStringSet"){
            X<- X[[1]]
        }
        results <- cbind(seqnames= 1,ExtractSubSequence(x=X,k=k,seq.size = seq.size))
    }
    X <- Biostrings::DNAStringSet(as.vector(results$seq))
    message("Check sequences composition...")
    resFreq <- Biostrings::letterFrequency(X,"N",as.prob = T)
    testNFreq <- as.vector(resFreq>0.1)
    if(any(testNFreq)){
        message(paste0("Warning: Some of your sequences have a N frequency > 0.1 and will be removed.\nDeepG4 has difficulty to handle sequences with a N rate > 10%"))

        X <- X[!testNFreq]
        results <- results[!testNFreq,]
        if(length(X)<1){
            stop("Not enough sequences to continue ...",
                 call. = FALSE)
        }
    }
    predictions <- DeepG4(X = X)
    results$score <- predictions[,1]
    results <- results[results$score>treshold,]
    if(nrow(results)== 0){
        stop(paste0("No sequences with a score <",treshold),
             call. = FALSE)
    }
    return(results)
}

#' Scanning of potential active G4 on a sequence or a list of sequences using a sliding window of size k.
#'
#' @param X An object of class GRanges.
#' @param X.ATAC An object of class GRanges containing DNA accessibility.
#' @param k size of the sliding windows.
#' @param treshold  numeric value who define the treshold to use to consider a sequence asc ontaining an active G4.
#' @param threads  numeric value who define the number of threads used in DeepG4Scan (Generate sub sequences)
#' @param lower.case a boolean. Set to \code{TRUE} if elements of X are in lower case (default to FALSE).
#' @param GENOME a BSgenome object containing the DNA sequence of genome of interest.
#' @param use.bg a boolean. Set to \code{TRUE} you want to normalize the accessibility using a windows background of windows_bg.
#' @param windows_bg numeric value who define the windows use to get background signal.
#' @param treshold_bg numeric value who set the treshold signal/background.
#'
#' @return a data.frame with the position of potential active G4 across input sequences.
#' @examples
.callDeepG4ScanATAC <- function(X=NULL,X.ATAC=NULL,k=20,treshold = 0.5,threads = 1,lower.case=F,GENOME = NULL,use.bg=T,windows_bg=5000,treshold_bg = 2){
    seq.size = 201
    # Packages check
    if (!requireNamespace("keras", quietly = TRUE)) {
        stop("Package \"keras\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (!requireNamespace("Biostrings", quietly = TRUE)) {
        stop("Package \"Biostrings\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("Package \"GenomicRanges\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (is.null(GENOME)) {
        stop("GENOME must be provided (see ?DeepG4Scan for accepted formats).",
             call. = FALSE)
    }
    if (class(GENOME) != "BSgenome") {
        stop("GENOME must be a BSgenome object.",
             call. = FALSE)
    }
    if (is.null(X.ATAC)) {
        stop("X.ATAC must be provided (see ?DeepG4Scan for accepted formats).",
             call. = FALSE)
    }
    if(class(X.ATAC) == "SimpleRleList"){
        X.ATAC <- as(X.ATAC,"GRanges")
    }else if(class(X.ATAC) == "GRanges"){
        if(!"score" %in% colnames(IRanges::values(X.ATAC))){
            stop("X.ATAC must be a GRanges object with a metadata column 'score'",
                 call. = FALSE)
        }
        if(class(X.ATAC$score) != "numeric"){
            stop("score column must be a numeric vector.",
                 call. = FALSE)
        }
    }else{
        stop("X.ATAC must be a character (path for a file in BedGraph/BigWig format), or a SimpleRleList/GRanges class",
             call. = FALSE)
    }
    seq.X <- unique(as.character(seqnames(X)))
    seq.X.atac <- unique(as.character(seqnames(X.ATAC)))
    X.in.atac <- !seq.X %in%seq.X.atac
    if(any(X.in.atac)){
        message(paste0("Warning: Some of your sequences have chromosome who's missing in the Accessibility file : ",paste(seq.X[X.in.atac],collapse=" "),", and will be croped."))
        X <- X[!as.character(seqnames(X)) %in% seq.X[X.in.atac]]
    }
    binbed <- rtracklayer::import.bed(system.file("extdata", "random_region_for_scaling_min_max.bed", package = "DeepG4"))
    message("Normalise accessibility ...")
    X.ATAC <- NormBW(X.ATAC,binbed)
    message(message("Extract sub-sequences / sub-accessibility using ",threads," threads..."))
    results <- ExtractSubSequencefromBed(x=X,x.atac=X.ATAC,GENOME=GENOME,nb.threads=threads,use.bg=T)
    message(paste0("From ",length(X)," positions to",nrow(results),"."))
    X <- Biostrings::DNAStringSet(as.vector(results$seq))
    x.atac <- results$score
    message("Check sequences composition...")
    resFreq <- Biostrings::letterFrequency(X,"N",as.prob = T)
    testNFreq <- as.vector(resFreq>0.1)
    if(any(testNFreq)){
        message(paste0("Warning: Some of your sequences have a N frequency > 0.1 and will be removed.\nDeepG4 has difficulty to handle sequences with a N rate > 10%"))
        X <- X[!testNFreq]
        x.atac <- x.atac[!testNFreq]
        if(length(X)<1){
            stop("Not enough sequences to continue ...",
                 call. = FALSE)
        }
    }
    ## Predictions using DeepG4
    message("prediction using DeepG4 main function...")
    predictions <- DeepG4(X = X,X.atac = x.atac,lower.case = lower.case)
    results <- results[!testNFreq,]
    results$score <- predictions[,1]
    results <- results[results$score>treshold,]
    if(nrow(results)== 0){
        stop(paste0("No sequences with a score <",treshold),
             call. = FALSE)
    }
    return(results)
}


#' Standard method to call
#'
#' @param X One of character, list, DNAString, DNAStringSet, DNAStringSetList, GRanges
#' @param ...
#'
#' @return
#' @export
#'
#' @details
#' This function is a method who launch ?callDeepG4Scan or .callDeepG4ScanATAC based on the arguments who are passed on. See ?.callDeepG4ScanATAC or ?.callDeepG4Scan.

setGeneric("DeepG4Scan", function(X,...) standardGeneric("DeepG4Scan"))

setMethod("DeepG4Scan",c(X="character"), .callDeepG4Scan)
setMethod("DeepG4Scan",c(X="list"), .callDeepG4Scan)
setMethod("DeepG4Scan",c(X="DNAString"), .callDeepG4Scan)
setMethod("DeepG4Scan",c(X="DNAStringSet"), .callDeepG4Scan)
setMethod("DeepG4Scan",c(X="DNAStringSetList"), .callDeepG4Scan)
setMethod("DeepG4Scan",c(X="GRanges"), .callDeepG4ScanATAC)
