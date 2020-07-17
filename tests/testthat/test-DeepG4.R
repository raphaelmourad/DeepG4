require(Biostrings)
library(rsample)
test_sequences <- Biostrings::readDNAStringSet(system.file("extdata", "promoters_seq_example.fa", package = "DeepG4"))
res.pos <- Biostrings::readDNAStringSet(system.file("extdata", "test_G4_data.fa", package = "DeepG4"))
shuffle <- function(dna) {
    # Convert to a vector of single bases
    charvec <- strsplit(as.character(dna),"")[[1]]
    # Shuffle the vector
    shuffled_charvec <- sample(charvec)
    # Convert back to a DNA string
    DNAString( paste(shuffled_charvec, collapse="") )
}

#Test errors

test_that("Test with NULL X", {

    expect_error(DeepG4())

})

test_that("Test with numerical X", {

    expect_error(DeepG4Scan(X=c(1,2,3)))

})

test_that("Test with list of numerical X", {

    expect_error(DeepG4Scan(X=list(c(1,2,3),c(1,2,3))))

})

test_that("One sequence with N>0.1", {
    myseq <- Biostrings::subseq(test_sequences[[1]],start = 1,width=201)
    samplesize <- round(Biostrings::nchar(myseq)*0.2)
    myseq <- replaceLetterAt(myseq,sample(1:Biostrings::nchar(myseq),samplesize),paste0(rep("N",samplesize),collapse=""))
    expect_error(DeepG4(X = myseq))
})


sequences.pos <- readDNAStringSet(system.file("extdata", "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b.Fa", package = "DeepG4"))
sequences.ctrl <- readDNAStringSet(system.file("extdata", "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM.Fa", package = "DeepG4"))
sequences <- c(sequences.pos,sequences.ctrl)
# Generate classes
Y <- c(rep(1,length(sequences.pos)),rep(0,length(sequences.ctrl)))
smp_size <- floor(0.70 * length(sequences))
train_ind <- sample(seq_len(length(sequences)), size = smp_size)
x.train <- sequences[train_ind]
x.test <- sequences[-train_ind]
y.train <- Y[train_ind]
y.test <- Y[-train_ind]

#Test with evaluation
test_that("Test with evaluation", {
    training <- DeepG4(x.train,y.train)
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

