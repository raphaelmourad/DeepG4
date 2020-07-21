require(Biostrings)
res.pos <- Biostrings::readDNAStringSet(system.file("extdata", "test_G4_data.fa", package = "DeepG4"))

#Test errors

test_that("Test with NULL X", {

    expect_error(ExtractMotifFromModel())

})

test_that("Test with numerical X", {

    expect_error(ExtractMotifFromModel(X=c(1,2,3)))

})

test_that("Test with list of numerical X", {

    expect_error(ExtractMotifFromModel(X=list(c(1,2,3),c(1,2,3))))

})

test_that("Test with one little sequence", {
    expect_is(ExtractMotifFromModel(DNAStringSet("AAAA")),"list")
})


test_that("Test with one little sequence", {
    expect_is(ExtractMotifFromModel(DNAStringSet(paste0(rep("A",20),collapse=""))),"list")
})

test_that("One sequence with N>0.1", {
    myseq <- Biostrings::subseq(res.pos[[1]],start = 1,width=201)
    samplesize <- round(Biostrings::nchar(myseq)*0.2)
    myseq <- replaceLetterAt(myseq,sample(1:Biostrings::nchar(myseq),samplesize),paste0(rep("N",samplesize),collapse=""))
    expect_error(ExtractMotifFromModel(X = myseq))
})


test_that("Test ExtractMotifFromModel with one seq", {
    res <- ExtractMotifFromModel(res.pos[[1]])
    expect_is(res, "list")
})

test_that("Two sequence with character object", {
    res <- ExtractMotifFromModel(X = as.character(res.pos[1:2]))
    expect_is(res, "list")
})

test_that("one sequence list with DNAString object", {
    res <- ExtractMotifFromModel(X = list(res.pos[[1]]))
    expect_is(res, "list")
})

test_that("Two sequence with a list object", {
    res <- ExtractMotifFromModel(X = lapply(res.pos[1:2],as.character))
    expect_is(res, "list")
})

test_that("Two sequence with DNAStringSetList object", {
    res <- ExtractMotifFromModel(X = DNAStringSetList(res.pos[1],res.pos[2]))
    expect_is(res, "list")
})

test_that("one sequence with seq.size > 201", {
    res <- ExtractMotifFromModel(X = paste0(res.pos[1],res.pos[1]))
    expect_is(res, "list")
})

test_that("one sequence with seq.size > 201", {
    res <- ExtractMotifFromModel(X = paste0(res.pos[1],res.pos[1]),top_kernel = 1000)
    expect_is(res, "list")
})
