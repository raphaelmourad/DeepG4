#' ExtractMotifFromModel allow the user to extract motifs who are activated by the sequences in our model.
#'
#' @param X An object of class character,list or DNAStringSet/DNAStringSetList with DNA sequences.
#' @param Y a numeric vector of 1 and 0 values (default to NULL).
#' @param lower.case boolean. Set to \code{TRUE} if elements of X are in lower case (default to FALSE).
#' @param top_kernel integer which indicate the top best motifs activated by X.
#' @details A kernel is a weight matrix used in convolution layer in DNN (Deep Neural Network).
#' With one-hot DNA, a kernel represented a weighted DNA motif, but can't be easily interpreted.
#' With this function, we detect the location on the response map where the kernel is activated for all input sequences and
#'  built the corresponding position count matrix.
#' @return A list of Position Count Matrix.
#' @export
ExtractMotifFromModel <- function(X = NULL,Y=NULL,lower.case=F,top_kernel = 20){
    seq.size <- 201
    kernel_size <- 20
    tabv = c("N"=5,"T"=4,"G"=3,"C"=2,"A"=1)
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
    }else if(class(X)[[1]] =="DNAString"){
        X <- Biostrings::DNAStringSet(X)
    }
    ## Check sequences sizes
    message("Check sequences sizes...")
    seqsizes <- Biostrings::nchar(X)
    if(max(seqsizes) > seq.size){
        message(paste0("Warning: Some of your sequences are >",seq.size,", and will be croped."))
        X[(Biostrings::nchar(X)>seq.size & !Biostrings::nchar(X)<seq.size)] <- Biostrings::subseq(X[(Biostrings::nchar(X)>seq.size & !Biostrings::nchar(X)<seq.size)],start = 1,end = seq.size)
    }
    ## Check DNA composition
    message("Check sequences composition...")
    resFreq <- Biostrings::letterFrequency(X,"N",as.prob = T)
    testNFreq <- as.vector(resFreq>0.1)
    if(any(testNFreq)){
        message(paste0("Warning: Some of your sequences have a N frequency > 0.1 and will be removed.\nDeepG4 has difficulty to handle sequences with a N rate > 10%"))
        X <- X[!testNFreq]
        if(length(X)<1){
            stop("Not enough sequences to continue ...",
                 call. = FALSE)
        }
    }
    ## One-Hot conversion
    message("One-Hot Conversion...")
    if(length(seqsizes) == 1) {
        X_oh <- DNAToNumerical(X,tabv = tabv,lower.case=lower.case,seq.size = seq.size)
    }else{
        ## Have to apply One-Hot independently because seq sizes are differents
        X_by_size <- lapply(unique(Biostrings::nchar(X)),function(onesize){
            DNAToNumerical(X[Biostrings::nchar(X)==onesize],tabv = tabv,lower.case=lower.case,seq.size = seq.size)
        })
        X_oh <- array(unlist(X_by_size), dim = c(length(X),seq.size,length(tabv)))
    }
    model <-  system.file("extdata", "model.hdf5", package = "DeepG4")
    #Load model with keras (tensorflow must be installed as well)
    model <- keras::load_model_hdf5(model)
    weights <- keras::get_weights(object = model)[[1]]
    Convolution <- keras::keras_model(inputs = model$input,
                               outputs = keras::get_layer(model, index = 2)$output)
    res <- stats::predict(Convolution,X_oh)
    kernels_information <- colSums(apply(res,c(1,3),max))
    nb_of_kernels <- length(kernels_information)
    if(top_kernel > nb_of_kernels){
        message(paste0("Warning: top_kernel > ",nb_of_kernels," (number of kernels)."))
        top_kernel <- nb_of_kernels
    }
    kernels_information <- order(kernels_information,decreasing = T)[1:top_kernel]
    if(max(res[,,kernels_information[1]]) == 0){
        stop(paste0("Error: No kernels activated, please try to select more sequences or larger (<=",seq.size,"bp)"),
             call. = FALSE)
    }
    PCM.list <- lapply(kernels_information,function(kernel){
        #Response map
        response_map <- res[,,kernel]
        #Make pwm from activated sequences
        if(is.null(dim(response_map))){
            max_values <- max(response_map)
            pos_in_seq <- which.max(response_map)
            end_in_seq <- ifelse(pos_in_seq+(kernel_size-1)<Biostrings::nchar(X),pos_in_seq+(kernel_size-1),Biostrings::nchar(X))
            DNA_seq <- Biostrings::subseq(X,start=pos_in_seq,end=end_in_seq)
        }else{
            max_values <- apply(response_map,1,max)
            indexes_no_0 <- max_values!=0
            DNA_seq <- X[indexes_no_0]
            if(length(DNA_seq)!=1){
                pos_in_seq <- apply(response_map[indexes_no_0,],1,which.max)
                end_in_seq <- pos_in_seq+(kernel_size-1)
                end_in_seq[end_in_seq>Biostrings::nchar(DNA_seq)] <- Biostrings::nchar(DNA_seq)


            }else{
                pos_in_seq <- which.max(response_map[indexes_no_0,])
                end_in_seq <- ifelse(pos_in_seq+(kernel_size-1)<Biostrings::nchar(DNA_seq),pos_in_seq+(kernel_size-1),Biostrings::nchar(DNA_seq))
            }
            DNA_seq <- Biostrings::subseq(DNA_seq,start=pos_in_seq,end=end_in_seq)

        }






        t(Reduce("+",
                 lapply(DNA_seq,function(x){Biostrings::letterFrequencyInSlidingView(x,letters = c("A","C","G","T"),
                                                                                     as.prob = F,view.width = 1)})
        ))
    })
    return(PCM.list)
}
