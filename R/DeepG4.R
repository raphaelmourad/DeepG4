DeepG4 <- function(sequences = NULL,model = NULL,tabv = c("N"=5,"T"=4,"G"=3,"C"=2,"A"=1),lower.case=F){
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
    if(!class(sequences) %in%c("DNAStringSet","DNAStringSetList")){
        if(class(sequences) == "character"){
            sequences <- Biostrings::DNAStringSet(sequences)
        }else if(class(sequences) == "list"){
            sequences <- unlist(Biostrings::DNAStringSetList(sequences))
        }else{
            stop("sequences must be a character, a list or a DNAStringSet/DNAStringSetList",
                 call. = FALSE)
        }
    }else if(class(sequences) =="DNAStringSetList"){
        sequences <- unlist(Biostrings::DNAStringSetList(sequences))
    }else{
        stop("sequences must be a character, a list or a DNAStringSet/DNAStringSetList",
             call. = FALSE)
    }
    ## Check sequences sizes
    Biostrings::width(sequences)
    ## One-Hot conversion
    sequences <- DNAToNumerical(sequences,tabv = tabv,lower.case=lower.case)
    # Try to load our saved model or custom model if !is.null(model)
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


}
