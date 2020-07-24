require(Biostrings)
test_sequences <- Biostrings::readDNAStringSet(system.file("extdata", "promoters_seq_example.fa", package = "DeepG4"))
res.pos <- Biostrings::readDNAStringSet(system.file("extdata", "test_G4_data.fa", package = "DeepG4"))

#Test errors

test_that("Test with NULL X", {

    expect_error(DeepG4())

})

test_that("Test with numerical X", {

    expect_error(DeepG4(X=c(1,2,3)))

})

test_that("Test with list of numerical X", {

    expect_error(DeepG4(X=list(c(1,2,3),c(1,2,3))))

})


test_that("One sequence bad Y (character)", {
    myseq <- Biostrings::subseq(test_sequences[[1]],start = 1,width=201)
    expect_error(DeepG4(X = myseq,Y="1"))
})

test_that("One sequence bad Y (with 1,2,3 values)", {
    myseq <- Biostrings::subseq(test_sequences[[1]],start = 1,width=201)
    expect_error(DeepG4(X = myseq,Y=c(1,2,3)))
})


test_that("One sequence bad Y (!= size)", {
    myseq <- Biostrings::subseq(test_sequences[[1]],start = 1,width=201)
    expect_error(DeepG4(X = myseq,Y=c(1,0)))
})



test_that("One sequence with N>0.1", {
    myseq <- Biostrings::subseq(test_sequences[[1]],start = 1,width=201)
    samplesize <- round(Biostrings::nchar(myseq)*0.2)
    myseq <- replaceLetterAt(myseq,sample(1:Biostrings::nchar(myseq),samplesize),paste0(rep("N",samplesize),collapse=""))
    expect_error(DeepG4(X = myseq))
})

test_that("Two sequence with N>=0.1", {
    myseq <- Biostrings::subseq(test_sequences[1:2],start = 1,width=201)
    samplesize <- round(Biostrings::nchar(myseq[1])*0.2)
    sampleMat <- matrix(FALSE,nrow = length(myseq),ncol = Biostrings::nchar(myseq[1]))
    for(i in 1:length(myseq)){
        sampleMat[i,sample(1:Biostrings::nchar(myseq[1]),samplesize)] <- TRUE
    }
    myseq <- replaceLetterAt(myseq,sampleMat,DNAStringSet(unlist(lapply(1:nrow(sampleMat),function(x){paste0(rep("N",samplesize),collapse="")}))))
    expect_error(DeepG4(X = myseq))
})

sequences.pos <- readDNAStringSet("/home/rochevin/Documents/PROJET_THESE/scriptsDeepG4/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b.Fa")[1:1000]
sequences.ctrl <- readDNAStringSet("/home/rochevin/Documents/PROJET_THESE/scriptsDeepG4/Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM.Fa")[1:1000]
sequences <- c(sequences.pos,sequences.ctrl)
# Generate classes
Y <- c(rep(1,length(sequences.pos)),rep(0,length(sequences.ctrl)))
smp_size <- floor(0.70 * length(sequences))
train_ind <- sample(seq_len(length(sequences)), size = smp_size)
x.train <- sequences[train_ind]
x.test <- sequences[-train_ind]
y.train <- Y[train_ind]
y.test <- Y[-train_ind]
#Test with training
test_that("Test with training no path", {
    training <- DeepG4(x.train,y.train,retrain=T)
    expect_is(training[[1]], "matrix")
    expect_is(training[[2]], "ggplot")
    expect_is(training[[3]], "ggplot")
    expect_is(training[[4]], "data.frame")
})
test_that("Test with training with path", {
    training <- DeepG4(x.train,y.train,retrain=T,retrain.path = "/home/rochevin/Documents/PROJET_THESE/DeepG4/model_retrained.h5")
    expect_is(training[[1]], "matrix")
    expect_is(training[[2]], "ggplot")
    expect_is(training[[3]], "ggplot")
    expect_is(training[[4]], "data.frame")
})
#Test with training
test_that("Test with training error Y", {
    expect_error(DeepG4(X=x.train,Y=c("cc1","cc2"),retrain=T))
    expect_error(DeepG4(X=x.train,Y=c(3,4),retrain=T))
    expect_error(DeepG4(X=x.train,retrain=T))
    expect_error(DeepG4(x.train,y.train,retrain.path="/home/rochevin/Documents/PROJET_THESE/DeepG4/model_retrained.h5",retrain=T,model=1))
})
#Test with evaluation
test_that("Test with evaluation", {
    training <- DeepG4(x.test,y.test,model = "/home/rochevin/Documents/PROJET_THESE/DeepG4/model_retrained.h5")
    expect_is(training[[1]], "matrix")
    expect_is(training[[2]], "ggplot")
    expect_is(training[[3]], "ggplot")
    expect_is(training[[4]], "data.frame")
})



test_that("One sequence with DNAStringSet object", {
    res <- DeepG4(X = test_sequences[1])
    expect_is(res, "matrix")
})

test_that("One sequence with DNAString object", {
    res <- DeepG4(X = test_sequences[[1]])
    expect_is(res, "matrix")
})

test_that("One sequence with character object", {
    res <- DeepG4(X = as.character(test_sequences[[1]]))
    expect_is(res, "matrix")
})

test_that("One sequence with a list object", {
    res <- DeepG4(X = list(as.character(test_sequences[[1]])))
    expect_is(res, "matrix")
})

test_that("One sequence with DNAStringSetList object", {
    res <- DeepG4(X = DNAStringSetList(test_sequences[1]))
    expect_is(res, "matrix")
})


test_that("One sequence with DNAStringSetList object", {
    res <- DeepG4(X = DNAStringSet(tolower(test_sequences[1])),lower.case = T)
    expect_is(res, "matrix")
})

#Test with two sequences



test_that("Two sequence with DNAStringSet object", {
    res <- DeepG4(X = test_sequences[1:2])
    expect_is(res, "matrix")
})


test_that("Two sequence with character object", {
    res <- DeepG4(X = as.character(test_sequences[1:2]))
    expect_is(res, "matrix")
})

test_that("Two sequence with a list object", {
    res <- DeepG4(X = lapply(test_sequences[1:2],as.character))
    expect_is(res, "matrix")
})


test_that("Two sequence with DNAStringSetList object", {
    res <- DeepG4(X = DNAStringSetList(test_sequences[1],test_sequences[2]))
    expect_is(res, "matrix")
})

test_that("Two sequence with N<0.1", {
    myseq <- Biostrings::subseq(test_sequences[1:2],start = 1,width=201)
    samplesize <- round(Biostrings::nchar(myseq[1])*0.05)
    sampleMat <- matrix(FALSE,nrow = length(myseq),ncol = Biostrings::nchar(myseq[1]))
    for(i in 1:length(myseq)){
        sampleMat[i,sample(1:Biostrings::nchar(myseq[1]),samplesize)] <- TRUE
    }
    myseq <- replaceLetterAt(myseq,sampleMat,DNAStringSet(unlist(lapply(1:nrow(sampleMat),function(x){paste0(rep("N",samplesize),collapse="")}))))
    res <- DeepG4(X = myseq)
    expect_is(res, "matrix")
})


test_that("Two sequence with exactly 201bp", {
    myseq <- lapply(test_sequences[1:2],Biostrings::subseq,start = 1,width=201)
    res <- DeepG4(X = myseq,treshold=0)
    expect_is(res, "matrix")
})


test_that("No control level (full 0 or full 1)", {
    myseq <- Biostrings::subseq(test_sequences[1:2],start = 1,width=201)
    res <- DeepG4(X = myseq,Y=c(0,0))
    expect_is(res, "matrix")
})


test_that("estimate levels = 0", {
    myseq <- Biostrings::subseq(test_sequences[1:2],start = 1,width=201)
    res <- DeepG4(X = myseq,Y=c(0,1))
    expect_is(res, "list")
})

test_that("Return log of the odds instead of probability (set log_odds to TRUE)", {
    res <- DeepG4(X = test_sequences[1:2],log_odds=T)
    expect_is(res, "matrix")
})
