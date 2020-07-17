#' DeepG4 main function to predict a probability to form a G4, given a DNA sequence.
#'
#' @param X An object of class character,list or DNAStringSet/DNAStringSetList with DNA sequences.
#' @param Y a numeric vector of 1 and 0 values (default to NULL).
#' @param model a path to a keras model in hdf5 format (default to NULL). Don't change it unless you want to use our function with a custom model.
#' @param tabv a named vector of numeric values representing the DNA to numerical conversion. Don't change it unless you want to use our function with a custom model.
#' @param lower.case boolean. Set to \code{TRUE} if elements of X are in lower case (default to FALSE).
#' @param seq.size numeric value representing the sequence size accepted by our model. Don't change it unless you want to use our function with a custom model.
#' @param treshold numeric value who define the treshold to use to get confusion matrix (default to 0.5).
#' @param retrain boolean. Set to \code{TRUE} if you want to retrain with your own dataset. Need Y to be provided (default to FALSE).
#' @param retrain.path file where retrained model will be saved.
#' @details
#'  This function is a wrapper to help people to get a prediction given any DNA sequence(s) of type ACGTN with our DeepG4 model.
#'  You don't have to use it to get a DeepG4 prediction, if you're familar with keras and tensorflow, you can access our model in hdf5 package using \code{system.file("extdata", "model.hdf5", package = "DeepG4")}.
#'  In complement, \code{\link{DNAToNumerical}} can help you to get the one-hot conversion needed by our model as input.
#'  If your sequences > \code{seq.size}, they will be cropped and sequences < \code{seq.size}, will be filled with zero padding.
#' @return if \code{Y = NULL}, return DeepG4 prediction for each value of X.
#'     if \code{Y} is provided, return a list with list(prediction for each value of X,a ggplot2 object representing AUC,a ggplot2 object representing confusion matrix,some metrics)
#' @export
#'
#' @examples
#' library(Biostrings)
#' library(DeepG4)
#'
#' sequences <- system.file("extdata", "test_G4_data.fa", package = "DeepG4")
#' sequences <- readDNAStringSet(sequences)
#'
#' predictions <- DeepG4(sequences)
#' head(predictions)
DeepG4 <- function(X = NULL,Y=NULL,model = NULL,tabv = c("N"=5,"T"=4,"G"=3,"C"=2,"A"=1),lower.case=F,seq.size = 201,treshold = 0.5,retrain=FALSE,retrain.path=""){
    model.regular.size.accepted <- 201
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
        X <- DNAStringSet(X)
    }else if(!class(X)[[1]] =="DNAStringSet"){
        stop("X must be a character, a list or a DNAStringSet/DNAStringSetList",
             call. = FALSE)
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
        if(length(X)<1){
            stop("Not enough sequences to continue ...",
                 call. = FALSE)
        }
        X <- X[!testNFreq]
        if(length(X)<1){
            stop("Not enough sequences to continue ...",
                 call. = FALSE)
        }
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
    if(retrain){
        # IF RETRAIN = TRUE
        message("retrain == TRUE")
        message("Model will be retrain using user input...")
        # Check Y
        if(is.null(Y)){
            stop("Y must be set if you want retrain our model.",
                 call. = FALSE)
        }
        if(class(Y) != "numeric"){
            stop("Y must be a numeric vector of 1 and 0 values",
                 call. = FALSE)
        }
        if(FALSE %in% (unique(Y) %in% c(0,1))){
            stop("Y must be a numeric vector of 1 and 0 values",
                 call. = FALSE)
        }
        # Build the model
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
        model <- keras::from_config(keras::get_config(model))
        keras::compile(model,
                       optimizer = 'rmsprop',
                       loss = 'binary_crossentropy',
                       metrics = list('accuracy')
        )
        # Retrain the model
        history <- keras::fit(model,
                              X,
                              Y,
                              epochs = 20,
                              batch_size = 128,
                              validation_split = 0.2,
                              verbose= 1)
        res <- stats::predict(model,X)
        if(retrain.path == ""){
            retrain.path <- paste0("DeepG4_retrained_",Sys.Date(),".hdf5")
        }
        keras::save_model_hdf5(model,retrain.path)
    }else{
        #IF RETRAIN = FALSE
        # Try to load our saved model or custom model if !is.null(model)
        message("Loading model...")
        if(is.null(model)){
            model <-  system.file("extdata", "model.hdf5", package = "DeepG4")
            if(seq.size != model.regular.size.accepted){
                message("Please don't manually set seq.size unless you want to use a custom model")
                seq.size <- model.regular.size.accepted
            }
        }else{
            if(class(model) != "character"){
                stop("model must be a path to a keras model in hdf5 format",
                     call. = FALSE)
            }
        }
        #Load model with keras (tensorflow must be installed as well)
        model <- keras::load_model_hdf5(model)
        res <- stats::predict(model,X)
    }
    # If Y is provided, instead of returning prediction, return accuracy / AUC
    if(is.null(Y)){
        return(res)
    }else{
        if(class(Y) != "numeric"){
            stop("Y must be a numeric vector of 1 and 0 values",
                 call. = FALSE)
        }
        if(FALSE %in% (unique(Y) %in% c(0,1))){
            stop("Y must be a numeric vector of 1 and 0 values",
                 call. = FALSE)
        }
        # Compute accuracy and AUC
        if (!requireNamespace("ggplot2", quietly = TRUE)) {
            stop("Package \"ggplot2\" needed for this function to work. Please install it.",
                 call. = FALSE)
        }
        if (!requireNamespace("yardstick", quietly = TRUE)) {
            stop("Package \"yardstick\" needed for this function to work. Please install it.",
                 call. = FALSE)
        }
        prediction_table <- data.frame(
            truth = as.factor(Y),
            pred_prob = res[,1]
        )
        prediction_table$estimate <- as.factor(ifelse(prediction_table$pred_prob<treshold,0,1))
        if(length(levels(prediction_table$truth))==1){
            prediction_table$truth <- factor(prediction_table$truth,levels = c(1,0))
        }
        #Plot AUC
        plot_ROC <- ggplot2::autoplot(yardstick::roc_curve(prediction_table,truth,pred_prob))
        #Get metrics
        table_metrics <- yardstick::metrics(prediction_table,truth,estimate)
        #Plot confusion matrix
        confusion_matrix <- as.data.frame(yardstick::conf_mat(prediction_table,truth, estimate)[[1]])

        confusion_matrix <- ggplot2::ggplot(confusion_matrix,ggplot2::aes(Prediction, Truth, fill = Freq)) +
            ggplot2::geom_tile(show.legend = FALSE) +
            ggplot2::scale_fill_viridis_c() +
            ggplot2::geom_text(ggplot2::aes(label = Freq), color = "white", alpha = 1, size = 8) +
            ggplot2::labs(
                title = "Confusion matrix"
            ) + ggplot2::theme_minimal(base_size=18)
        return(list(res,plot_ROC,confusion_matrix,table_metrics))
    }


}
