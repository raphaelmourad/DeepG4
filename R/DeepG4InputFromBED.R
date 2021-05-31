#' Title
#'
#' @param BED
#' @param ATAC
#' @param is.bw
#' @param GENOME
#'
#' @return
#' @export
#'
#' @examples
DeepG4InputFromBED <- function(BED = NULL,ATAC = NULL,is.bw = TRUE,GENOME = NULL){
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
        if(!"score" %in% colnames(IRanges::values(ATAC))){
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
    binbed <- rtracklayer::import.bed(system.file("extdata", "random_region_for_scaling_min_max.bed", package = "DeepG4"))
    ATAC <- NormBW(ATAC,binbed)
    X.ATAC <- getScoreBW(ATAC,BED)
    X <- Biostrings::getSeq(GENOME,BED)
    return(list(X,as.vector(X.ATAC)))
}
