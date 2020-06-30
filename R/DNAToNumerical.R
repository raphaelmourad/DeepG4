DNAToNumerical <- function(x,tabv = c("N"=5,"T"=4,"G"=3,"C"=2,"A"=1),lower.case=F){
    if(lower.case){
        names(tabv) <- tolower(tabv)
    }
    x <- as.matrix(x)
    listMat=list()
    for(i in tabv){
        mat=matrix(0,nrow(x),ncol(x))
        mat[x==names(i)]=1
        listMat[[i]]=mat
    }
    arrayout=array(unlist(listMat), dim = c(nrow(listMat[[1]]), ncol(listMat[[1]]), length(listMat)))
    return(arrayout)

}
