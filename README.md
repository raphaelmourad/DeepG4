
<!-- README.md is generated from README.Rmd. Please edit that file -->

![logo](logo.svg)

<!-- ### [__DeepG4__: A deep learning approach to predict active G-quadruplexes](https://www.biorxiv.org/content/early/2020/07/23/2020.07.22.215699) -->

### **DeepG4**: A deep learning approach to predict cell-type specific active G-quadruplex regions

*Vincent Rocher, Matthieu Genais, Elissar Nassereddine and Raphael
Mourad*

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/morphos30/DeepG4/branch/master/graph/badge.svg)](https://codecov.io/gh/morphos30/DeepG4?branch=master)
<!-- badges: end -->

**DeepG4** is a deep learning model developed to predict a score of DNA
sequences to form active G-Guadruplexes (found both in vitro and in
vivo) using **DNA sequences** and **DNA accessibility**. **DeepG4** is
built in keras+tensorflow and is wrapped in an R package.

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

### With accessibility

Given small regions (bed) and an accessibility file (coverage file from
ATAC-seq/DNAse-seq/MNase-seq), you can predict active G4 regions in a
**specific cell type**:

``` r
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(DeepG4)

BED <- system.file("extdata", "test_G4_data.bed", package = "DeepG4")
BED <- import.bed(BED)
ATAC <- system.file("extdata", "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Accessibility.bw", package = "DeepG4")
ATAC <- import.bw(ATAC)

Input_DeepG4 <- DeepG4InputFromBED(BED=BED,ATAC = ATAC,GENOME=BSgenome.Hsapiens.UCSC.hg19)

Input_DeepG4
```

    [[1]]
      A DNAStringSet instance of length 100
          width seq
      [1]   201 GTTCGGGCCTCGGTCGCGCCGCCGGGTCTTGCAGACGCGAATGTAAACAGAAACA...TGACTCCTGGAGCGACCTTCACGAGGGAAAGCGCGCCCCCCGGCACCCACCCCT
      [2]   201 TTTCTATAGTTTTCTTTTGTTTCTACCTCATGACTAGATGATTCACTGCTTGAAC...GTCAAATCTGTCCATCTTCACTGCCACCCTTCAGTACCAAATGACCAGTCTCTT
      [3]   201 GCTTAAAAGCCTGTAAGAAAGATATAATTTGATAGAACTGGCTAGGATTTGTCAG...CGTCAGGGAGGGGGTGGGGCCTCCACGTGGGAGATCTTGCCTGGAGGTGGTGGA
      [4]   201 TCCCACACCCGGTAGATGTAAGGGAAAAACTGCATTACCCAGAAGGCACTGCCCC...GTGTGACGTCATCTCCGTGGGCCGGTTTGGCCCTGAAACAGTGTGGGGCCTAGA
      [5]   201 AGTAGCTACAGAGTTCCTGCTCCAGCAACCAGGAGCCTTGAGGCAGCACAAGGAC...ACCACAATGTCTGCCAAGAAAGAGGATGAGTCACCAAGACCCACAGGAAAGAGG
      ...   ... ...
     [96]   201 CACATGCCTTCCTTGGGGACGTGTTCACACATGTGGCCCTAGCTGTGAGAGACAG...CATCTCAGAACAGCTGAGCTGGAAGTGGGTGAATAATAATAATAATAATAATAA
     [97]   201 TGGTGGTCTTTCTCTACCGGGCCTGGTAGCCAAAGACAAAGGTCATAATCACTTG...CTATGTACTCTTCAAAGTGCCACCTCCTGGCTGCAAGCCAACCAACACAAAACC
     [98]   201 TGACCGTAGACCTCGTGCACTTCTGCTGCGGTCGGGGCCGGAGTCTGGGCTGGAG...GCGATCCAGAGCCAAGCGCCCCGCCCCTGCCCGGGCGCGCTCCCTCCTTAGCCC
     [99]   201 TTAACGTCATCAGTCGGGAGGACGACAGCTACGCACGCGCGGGGCACCTCCTCTG...GCCACGGTGGAGGCAGCGGCGAGAGGGGGCGGGGACAAGGAGAGGGCACGCACG
    [100]   201 GTGTCCGGGTGAGAGACCTGGAGGTGGGGCCTAGGTGTCTACCCGGCCAGGTGCG...TAAGGCTCGGGGCCAGTCGTCGTCCATTCCTTCCTAACACCTCCCTATCCTCCC

    [[2]]
      [1] 0.000000000 0.016287416 0.033261447 0.069375103 0.018520650 0.010934717 0.036308476 0.315843234 0.037658374
     [10] 0.045887551 0.037320211 0.042853401 0.068908093 0.071774485 0.084947561 0.027456211 0.033915868 0.006912598
     [19] 0.012604675 0.051405275 0.093813195 0.019288668 0.051228826 0.019520666 0.048686840 0.050116329 0.045801884
     [28] 0.033079207 0.035834917 0.056326946 0.096531489 0.064706374 0.026422647 0.016979087 0.008512502 0.021891554
     [37] 0.016688682 0.109472225 0.047901838 0.066676075 0.052591085 0.017467983 0.035541899 0.060001992 0.028878783
     [46] 0.056284886 0.045126048 0.052469122 0.101620595 0.047741155 0.036925371 0.021645371 0.044472962 0.012457179
     [55] 0.020373459 0.109529076 0.039006694 0.047824384 0.028752257 0.015437852 0.069926660 0.022213134 0.019726120
     [64] 0.044609840 0.028773493 0.008077349 0.042587371 0.016502886 0.035757895 0.015023933 0.024181422 0.057516040
     [73] 0.027492004 0.030316917 0.049878433 0.020105394 0.025934350 0.023845766 0.032338052 0.048007935 0.136436151
     [82] 0.060423998 0.034617445 0.051958662 0.064664156 0.034518694 0.020277026 0.042060108 0.055335700 0.051632313
     [91] 0.066588875 0.030586623 0.043823259 0.034947155 0.082091662 0.008496193 0.034567766 0.055516400 0.062191534
    [100] 0.049011882

Then predict using both **DNA** and **Accessibility** :

``` r
predictions <- DeepG4(X=Input_DeepG4[[1]],X.atac = Input_DeepG4[[2]])
head(predictions)
```

              [,1]
    [1,] 0.2159184
    [2,] 0.8819393
    [3,] 0.9991976
    [4,] 0.9999995
    [5,] 0.9740031
    [6,] 0.2259631

### Without accessbility

You still can predict active G4 regions using only **DNA** sequences :

``` r
library(rtracklayer)
library(DeepG4)

sequences <- readDNAStringSet(system.file("extdata", "test_G4_data.fa", package = "DeepG4"))
predictions <- DeepG4(X=sequences)
head(predictions)
```

              [,1]
    [1,] 0.9478214
    [2,] 0.5868858
    [3,] 0.9660227
    [4,] 0.9093548
    [5,] 0.9119551
    [6,] 0.2471965

<!-- ## Advanced usage of DeepG4 -->
<!-- If you have a large sequence (>201bp up to several Mbp), you can scan the sequence  and predict the positions of active G4s within the sequence. -->
<!-- ``` r -->
<!-- library(Biostrings) -->
<!-- library(DeepG4) -->
<!-- sequences <- readDNAStringSet(system.file("extdata", "promoters_seq_example.fa", package = "DeepG4")) -->
<!-- res <- DeepG4Scan(X = sequences,k=20,treshold=0.5) -->
<!-- ``` -->
<!-- DeepG4Scan function scans each input sequence with a step of  `k=20` and outputs for each input sequence the G4 positions (+/- 100bp) and the corresponding DeepG4 probabilities (>= treshold). -->
<!-- ``` r -->
<!-- library(dplyr) -->
<!-- res %>% dplyr::select(-seq) %>% group_by(seqnames) %>% dplyr::slice(1:2) %>%  head -->
<!-- ``` -->
<!-- ``` -->
<!-- # A tibble: 6 x 5 -->
<!-- # Groups:   seqnames [3] -->
<!--   seqnames start   end width score -->
<!--      <int> <int> <int> <int> <dbl> -->
<!-- 1        1  1241  1441   201 0.670 -->
<!-- 2        1  1261  1461   201 0.659 -->
<!-- 3        2  1481  1681   201 0.648 -->
<!-- 4        2  1521  1721   201 0.517 -->
<!-- 5        3  2161  2361   201 0.723 -->
<!-- 6        3  2181  2381   201 0.998 -->
<!-- ``` -->

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
library(cowplot)
sequences <- readDNAStringSet(system.file("extdata", "test_G4_data.fa", package = "DeepG4"))
res <- ExtractMotifFromModel(sequences,top_kernel=4)
p.pcm <- lapply(res,function(x){ggseqlogo(as.matrix(x)) + ggplot2::theme_classic(base_size=14)})
print(plot_grid(plotlist = p.pcm,ncol=2))
```

![](best_pcm_from_kernel.svg)

## Using DeepG4 with a new active G4 dataset

If you want to use our model architecture, but retrain with your own
dataset, you can do it by running our function `DeepG4` with
`retrain = TRUE`

``` r
library(Biostrings)
library(DeepG4)
library(rsample)
library(BSgenome.Hsapiens.UCSC.hg19)


ATAC <- system.file("extdata", "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Accessibility.bw", package = "DeepG4")
ATAC <- import.bw(ATAC)
# Read positive and segative set of sequences 
bed.pos <- import.bed(system.file("extdata", "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b.bed", package = "DeepG4"))
bed.neg <- import.bed(system.file("extdata", "Peaks_BG4_G4seq_HaCaT_GSE76688_hg19_201b_Ctrl_gkmSVM.bed", package = "DeepG4"))

# Generate classes
Y <- c(rep(1,length(bed.pos)),rep(0,length(bed.neg)))
BED <- c(bed.pos,bed.neg)
Input_DeepG4 <- DeepG4InputFromBED(BED=BED,ATAC = ATAC,GENOME=BSgenome.Hsapiens.UCSC.hg19)
```

``` r
training <- DeepG4(X=Input_DeepG4[[1]],X.atac=Input_DeepG4[[2]],Y,retrain=TRUE,retrain.path = "DeepG4_retrained.hdf5")
```

You can now take a look on the results :

``` r
library(cowplot)
p_res_train <- cowplot::plot_grid(plotlist = training[2:3])
print(p_res_train)
```

![](p_res_train.svg)

``` r
training[[4]]
```

    # A tibble: 4 x 3
      .metric     .estimator .estimate
      <chr>       <chr>          <dbl>
    1 accuracy    binary        0.989 
    2 kap         binary        0.978 
    3 mn_log_loss binary        0.0429
    4 roc_auc     binary        0.999 
