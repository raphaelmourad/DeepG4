#' Scanning of potential active G4 on a sequence or a list of sequences using a sliding window of size k.
#'
#' @param X An object of class character,list or DNAStringSet/DNAStringSetList with DNA sequences.
#' @param k size of the sliding windows.
#' @param treshold  numeric value who define the treshold to use to consider a sequence asc ontaining an active G4.
#' @param threads  numeric value who define the number of threads used in DeepG4Scan (Generate sub sequences)
#'
#' @return a data.frame with the position of potential active G4 across input sequences.
#' @export
#'
#' @examples
DeepG4Scan <- function(X=NULL,k=20,treshold = 0.5,threads = 1){
    seq.size = 201
    #Check if X is provided
    if (is.null(X)) {
        stop("X must be provided (see ?DeepG4 for accepted formats).",
             call. = FALSE)
    }
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


