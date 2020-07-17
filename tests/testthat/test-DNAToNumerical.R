require(Biostrings)
test_sequences <- Biostrings::readDNAStringSet(system.file("extdata", "promoters_seq_example.fa", package = "DeepG4"))
test_that("Test DNAToNumerical", {
    res <- DNAToNumerical(x = test_sequences[1])
    expect_is(res, "array")
})


test_that("Test DNAToNumerical tolower", {
    res <- DNAToNumerical(x = DNAStringSet(tolower(test_sequences[1])),lower.case = T)
    expect_is(res, "array")
})


test_that("Test DNAToNumerical seq size < 201", {
    myseq <- Biostrings::subseq(test_sequences[1:2],start = 1,width=100)
    res <- DNAToNumerical(x = myseq)
    expect_is(res, "array")
})
