
<!-- README.md is generated from README.Rmd. Please edit that file -->

![logo](logo.svg)

# **DeepG4**: A deep learning approach to predict active G-quadruplexes

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/morphos30/DeepG4/branch/master/graph/badge.svg)](https://codecov.io/gh/morphos30/DeepG4?branch=master)
<!-- badges: end -->

**DeepG4** is a deep learning model developed to predict the probability
of DNA sequences to form active G-Guadruplexes (found both in vitro and
in vivo). **DeepG4** is built in keras+tensorflow and is wrapped in an R
package.

## Requirements

DeepG4 was built with `Keras 2.3.1` and `tensorflow 2.1.0`, but it
should work with any version of theses libraries.

A very convenient way to install keras and tensorflow is using `R`. The
command line to install is from : <https://keras.rstudio.com/>.

``` r
install.packages("keras")
library(keras)
install_keras()
```

This will provide you with default CPU installations of Keras and
TensorFlow python packages (within a virtualenv) that can be used with
or without R.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("morphos30/DeepG4")
```

## Basic usage of DeepG4

If you have a small sequence (201bp or less), you can predict the
probability that the sequence forms an active G4.

``` r
library(Biostrings)
library(DeepG4)

sequences <- system.file("extdata", "test_G4_data.fa", package = "DeepG4")
sequences <- readDNAStringSet(sequences)

predictions <- DeepG4(sequences)
head(predictions)
```

``` 
          [,1]
[1,] 0.9998598
[2,] 0.9993761
[3,] 0.9539083
[4,] 0.9974855
[5,] 0.9908580
[6,] 0.9999917
```

## Advanced usage of DeepG4

If you have a large sequence (\>201bp up to several Mbp), you can scan
the sequence and predict the positions of active G4s within the
sequence.

``` r
library(Biostrings)
library(DeepG4)
sequences <- readDNAStringSet(system.file("extdata", "promoters_seq_example.fa", package = "DeepG4"))
res <- DeepG4Scan(X = sequences,k=20,treshold=0.5)
```

DeepG4Scan function scans each input sequence with a step of `k=20` and
outputs for each input sequence the G4 positions (+/- 100bp) and the
corresponding DeepG4 probabilities (\>= treshold).

``` r
library(dplyr)
res %>% dplyr::select(-seq) %>% group_by(seqnames) %>% dplyr::slice(1:2) %>%  head
```

    # A tibble: 6 x 5
    # Groups:   seqnames [3]
      seqnames start   end width score
         <int> <int> <int> <int> <dbl>
    1        1  1241  1441   201 0.670
    2        1  1261  1461   201 0.659
    3        2  1481  1681   201 0.648
    4        2  1521  1721   201 0.517
    5        3  2161  2361   201 0.723
    6        3  2181  2381   201 0.998

## SNP effect on g-quadruplex using DeepG4

Using our model, you can predict the potential effect of a SNP on active
G4 formation :

``` r
GetSeqFromSNPs <- function(my_granges,wsize = 201){
    SNP_pos <- (wsize - 1)/2 + 1 
    ## Compute Fasta
    SNps.seq.ref <- my_granges %>% anchor_center() %>% mutate(width = wsize) %>% getSeq(BSgenome.Hsapiens.UCSC.hg19.masked,.)
    ## Replace ref by alt
    sampleMat <- matrix(FALSE,nrow = length(SNps.seq.ref),ncol = nchar(SNps.seq.ref[1]))
    sampleMat[,SNP_pos] <- TRUE
    SNps.seq.alt <- replaceLetterAt(SNps.seq.ref, sampleMat, my_granges$alt)
    return(c(SNps.seq.ref,SNps.seq.alt))
}
require(GenomicRanges)
require(Biostrings)
require(dplyr)
require(plyranges)
require(BSgenome.Hsapiens.UCSC.hg19.masked)
SNPs <- GRanges(c("chr16:87350773","chr19:50093572"))
SNPs$name <- c("rs3748393","rs7249925")
SNPs$ref <- c("C","A")
SNPs$alt <- c("A","G")

SNPs_seq <- SNPs %>% GetSeqFromSNPs

DeepG4.score <- DeepG4(SNPs_seq)
SNPs$DeepG4_ref <- DeepG4.score[1:length(SNPs),]
SNPs$DeepG4_alt <- DeepG4.score[(length(SNPs)+1):nrow(DeepG4.score),]
SNPs <- SNPs %>% mutate(DeltaScore = DeepG4_alt-DeepG4_ref)
SNPs
```

    GRanges object with 2 ranges and 5 metadata columns:
          seqnames    ranges strand |         ref         alt        DeepG4_ref
             <Rle> <IRanges>  <Rle> | <character> <character>         <numeric>
      [1]    chr16  87350773      * |           C           A 0.840066134929657
      [2]    chr19  50093572      * |           A           G 0.126904487609863
                 DeepG4_alt         DeltaScore
                  <numeric>          <numeric>
      [1] 0.386543124914169 -0.453523010015488
      [2] 0.642100036144257  0.515195548534393
      -------
      seqinfo: 2 sequences from an unspecified genome; no seqlengths

## Scan DeepG4 DNA motifs from the input sequences

Using one-hot encoding of DNA, convolution kernels (first layer of
DeepG4) can be interpreted as weighted motifs, similar to position
weight matrices (PWMs) used for DNA motifs. The function
ExtractMotifFromModel detects DeepG4 DNA motifs found in the input
sequences.

``` r
library(Biostrings)
library(DeepG4)
library(ggseqlogo)
sequences <- readDNAStringSet(system.file("extdata", "test_G4_data.fa", package = "DeepG4"))
res <- ExtractMotifFromModel(sequences)
p.pcm <- lapply(res,function(x){ggseqlogo(as.matrix(x)) + ggplot2::theme_classic(base_size=8)})
print(cowplot::plot_grid(plotlist = p.pcm,ncol=4))
```

![](best_pcm_from_kernel.svg)
