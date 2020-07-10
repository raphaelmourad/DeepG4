#' Scanning of potential active G4 given a sequence or a list of sequences.
#'
#' @param X An object of class character,list or DNAStringSet/DNAStringSetList with DNA sequences.
#' @param seq.size numeric value representing the sequence size accepted by our model. Don't change it unless you want to use our function with a custom model.
#' @param treshold  numeric value who define the treshold to use to consider a sequence asc ontaining an active G4.
#' @param model a path to a keras model in hdf5 format (default to NULL). Don't change it unless you want to use our function with a custom model.
#'
#' @return a data.frame with the position of potential active G4 across input sequences.
#' @export
#'
#' @examples
DeepG4Scan <- function(X=NULL,seq.size = 201,treshold = 0.5,model = NULL){
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
            X <- unlist(Biostrings::DNAStringSetList(X))
        }else{
            stop("X must be a character, a list or a DNAStringSet/DNAStringSetList class",
                 call. = FALSE)
        }
    }else if(class(X)[[1]] =="DNAStringSetList"){
        X <- unlist(Biostrings::DNAStringSetList(X))
    }else{
        stop("X must be a character, a list or a DNAString/DNAStringSet/DNAStringSetList",
             call. = FALSE)
    }
    if(length(X)>1){
        results <- lapply(1:length(X),function(i){
            x <- X[i][[1]]
            start <- seq(1,length(x),seq.size)
            end <- start + (seq.size-1)
            end <- ifelse(end>length(x),length(x),end)
            Viewseq <- Biostrings::Views(x, start=start, end=end)
            sequences <- Biostrings::DNAStringSet(Viewseq)
            predictions <- DeepG4(X = sequences,seq.size = seq.size,model = model)

            results <- cbind(index = i,as.data.frame(IRanges::ranges(Viewseq)),seq=as.character(Viewseq),predictions = predictions[,1])
            results <- results[results$predictions>treshold,]
            return(results)
        })
        results <- do.call(rbind,results)
    }else{
        if(class(X)[[1]] =="DNAStringSet"){
            X<- X[[1]]
        }
        start <- seq(1,length(X),seq.size)
        end <- start + (seq.size-1)
        end <- ifelse(end>length(X),length(X),end)
        Viewseq <- Biostrings::Views(X, start=start, end=end)
        sequences <- Biostrings::DNAStringSet(Viewseq)
        predictions <- DeepG4(X = sequences,seq.size = seq.size,model = model)

        results <- cbind(as.data.frame(IRanges::ranges(Viewseq)),seq=as.character(Viewseq),predictions = predictions[,1])
        results <- results[results$predictions>treshold,]
    }
    return(results)
}


