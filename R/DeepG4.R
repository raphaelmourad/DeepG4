DeepG4 <- function(X = NULL,Y=NULL,model = NULL,tabv = c("N"=5,"T"=4,"G"=3,"C"=2,"A"=1),lower.case=F,seq.size = 201){
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
    # Check sequences and convert into one-hot
    ## Check model class and convert into DNAStringSet object if needed
    if(!class(X) %in%c("DNAStringSet","DNAStringSetList")){
        if(class(X) == "character"){
            X <- Biostrings::DNAStringSet(X)
        }else if(class(X) == "list"){
            X <- unlist(Biostrings::DNAStringSetList(X))
        }else{
            stop("X must be a character, a list or a DNAStringSet/DNAStringSetList class",
                 call. = FALSE)
        }
    }else if(class(X) =="DNAStringSetList"){
        X <- unlist(Biostrings::DNAStringSetList(X))
    }else{
        stop("X must be a character, a list or a DNAStringSet/DNAStringSetList",
             call. = FALSE)
    }
    ## Check sequences sizes
    message("Check sequences sizes...")
    seqsizes <- unique(Biostrings::nchar(X))
    if(max(seqsizes) > seq.size){
        message(paste0("Warning:: Some of your sequences are >",seq.size,", and will be croped."))
        X[(Biostrings::nchar(X)>seq.size & !Biostrings::nchar(X)<seq.size)] <- Biostrings::subseq(X[(Biostrings::nchar(X)>seq.size & !Biostrings::nchar(X)<seq.size)],start = 1,end = seq.size)
    }
    ## One-Hot conversion
    message("One-Hot Conversion...")
    if(length(seqsizes) == 1) {
        X <- DNAToNumerical(X,tabv = tabv,lower.case=lower.case,seq.size = seq.size)
    }else{
        ## Have to apply One-Hot independently because seq sizes are differents
        X_by_size <- lapply(unique(Biostrings::nchar(X)),function(onesize){
            DNAToNumerical(X[Biostrings::nchar(X)==onesize],tabv = tabv,lower.case=lower.case,seq.size = seq.size)
        })
        X <- array(unlist(X_by_size), dim = c(length(X),seq.size,length(tabv)))
    }
    # Try to load our saved model or custom model if !is.null(model)
    message("Loading model...")
    if(is.null(model)){
        model <-  system.file("extdata", "model.hdf5", package = "DeepG4")
    }else{
        if(class(model) != "character"){
            stop("model must be a path to a keras model in hdf5 format",
                 call. = FALSE)
        }
    }
    #Load model with keras (tensorflow must be installed as well)
    model <- keras::load_model_hdf5(model)
    res <- stats::predict(model,X)
    #If Y is provided, instead of returning prediction, return accuracy / AUC
    if(!is.null(Y)){
        if(class(Y) != "numeric"){
            stop("Y must be a numeric vector of 1 and 0 values",
                 call. = FALSE)
        }
        if(!unique(Y) %in% c(0,1)){
            stop("Y must be a numeric vector of 1 and 0 values",
                 call. = FALSE)
        }

    }

    return(res)
}
