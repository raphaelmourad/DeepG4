require(Biostrings)
test_sequences <- Biostrings::readDNAStringSet(system.file("extdata", "promoters_seq_example.fa", package = "DeepG4"))

#Test errors

test_that("Test with NULL X", {

    expect_error(DeepG4Scan())

})

test_that("Test with numerical X", {

    expect_error(DeepG4Scan(X=c(1,2,3)))

})

test_that("Test with list of numerical X", {

    expect_error(DeepG4Scan(X=list(c(1,2,3),c(1,2,3))))

})

test_that("One sequence with N>0.1", {
    myseq <- Biostrings::subseq(test_sequences[[1]],start = 1,width=1000)
    samplesize <- round(Biostrings::nchar(myseq)*0.2)
    myseq <- replaceLetterAt(myseq,sample(1:Biostrings::nchar(myseq),samplesize),paste0(rep("N",samplesize),collapse=""))
    expect_error(DeepG4Scan(X = myseq))
})


test_that("One sequence with score < treshold", {
    myseq <- Biostrings::subseq(test_sequences[[1]],start = 1,width=201)
    expect_error(DeepG4Scan(X = myseq))
})


#One sequence

test_that("One sequence with DNAStringSet object", {
    res <- DeepG4Scan(X = test_sequences[1])
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})

test_that("One sequence with DNAString object", {
    res <- DeepG4Scan(X = test_sequences[[1]])
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})

test_that("One sequence with character object", {
    res <- DeepG4Scan(X = as.character(test_sequences[[1]]))
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})

test_that("One sequence with a list object", {
    res <- DeepG4Scan(X = list(as.character(test_sequences[[1]])))
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})

test_that("One sequence with DNAStringSetList object", {
    res <- DeepG4Scan(X = DNAStringSetList(test_sequences[1]))
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})


test_that("One sequence with N<0.1", {
    myseq <- test_sequences[[1]]
    samplesize <- round(Biostrings::nchar(myseq)*0.05)
    myseq <- replaceLetterAt(myseq,sample(1:Biostrings::nchar(myseq),samplesize),paste0(rep("N",samplesize),collapse=""))
    res <- DeepG4Scan(X = myseq)
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})


test_that("One sequence with exactly 201bp", {
    myseq <- Biostrings::subseq(test_sequences[[1]],start = 1,width=201)
    res <- DeepG4Scan(X = myseq,treshold=0)
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})


#Two sequence

test_that("Two sequence with DNAStringSet object", {
    res <- DeepG4Scan(X = test_sequences[1:2])
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})


test_that("Two sequence with character object", {
    res <- DeepG4Scan(X = as.character(test_sequences[1:2]))
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})

test_that("Two sequence with a list object", {
    res <- DeepG4Scan(X = lapply(test_sequences[1:2],as.character))
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})


test_that("Two sequence with DNAStringSetList object", {
    res <- DeepG4Scan(X = DNAStringSetList(test_sequences[1],test_sequences[2]))
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})

test_that("Two sequence with N>=0.1", {
    myseq <- test_sequences[1:2]
    samplesize <- round(Biostrings::nchar(myseq[1])*0.1)
    sampleMat <- matrix(FALSE,nrow = length(myseq),ncol = Biostrings::nchar(myseq[1]))
    for(i in 1:length(myseq)){
        sampleMat[i,sample(1:Biostrings::nchar(myseq[1]),samplesize)] <- TRUE
    }
    myseq <- replaceLetterAt(myseq,sampleMat,DNAStringSet(unlist(lapply(1:nrow(sampleMat),function(x){paste0(rep("N",samplesize),collapse="")}))))
    res <- DeepG4Scan(X = myseq)
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})

test_that("Two sequence with N<0.1", {
    myseq <- test_sequences[1:2]
    samplesize <- round(Biostrings::nchar(myseq[1])*0.05)
    sampleMat <- matrix(FALSE,nrow = length(myseq),ncol = Biostrings::nchar(myseq[1]))
    for(i in 1:length(myseq)){
        sampleMat[i,sample(1:Biostrings::nchar(myseq[1]),samplesize)] <- TRUE
    }
    myseq <- replaceLetterAt(myseq,sampleMat,DNAStringSet(unlist(lapply(1:nrow(sampleMat),function(x){paste0(rep("N",samplesize),collapse="")}))))
    res <- DeepG4Scan(X = myseq)
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})


test_that("Two sequence with exactly 201bp", {
    myseq <- lapply(test_sequences[1:2],Biostrings::subseq,start = 1,width=201)
    res <- DeepG4Scan(X = myseq,treshold=0)
    expect_named(res,c("seqnames","start","end","width","seq","score"))
    expect_is(res, "data.frame")
})
