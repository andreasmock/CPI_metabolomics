CPI metabolomics data set
================
Andreas Mock, National Center for Tumor Diseases (NCT) Heidelberg

Load required R packages
========================

``` r
require(SummarizedExperiment)
```

Load data set
=============

``` r
load("CPI_metabolomics.RData")
```

Show object summary information:

``` r
CPI_metabolomics
```

    ## class: SummarizedExperiment 
    ## dim: 134 62 
    ## metadata(0):
    ## assays(1): ''
    ## rownames(134): C0 C2 ... SM C26:1 H1
    ## rowData names(8): id class ... fraction bond_type
    ## colnames: NULL
    ## colData names(3): X Gender Age

The metabolomic data is stored alongside the metabolite annotation and clinical annotation within a `SummarizedExperiment` class object. For more information regarding the see the corresponding [package tutorial](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html).

Access slots with data object
=============================

Access metabolomic data

``` r
head(assay(CPI_metabolomics))
```

    ##            [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
    ## C0    56.178169 42.680569 34.513216 28.848106 36.884648 39.041752
    ## C2     7.516536 11.703184  7.144435  7.711436  6.100314  3.178401
    ## C3     0.533766  0.415079  0.469999  0.452283  0.583664  0.371191
    ## C4     0.307557  0.268847  0.688515  0.532588  0.581600  0.199104
    ## C14:1  0.245055  0.210303  0.174569  0.154334  0.133054  0.164003
    ## C16    0.162439  0.219469  0.119258  0.109587  0.096276  0.057504
    ##            [,7]      [,8]      [,9]     [,10]     [,11]     [,12]
    ## C0    49.919588 39.767304 44.037454 43.538876 39.823094 32.306491
    ## C2     6.741930  4.015577  8.410512  3.975978  9.373824  9.314829
    ## C3     0.529669  0.393925  0.398677  0.384065  0.353552  0.387238
    ## C4     0.265489  0.126431  0.195943  0.201190  0.187512  0.959654
    ## C14:1  0.116527  0.152774  0.144088  0.166176  0.230703  0.195277
    ## C16    0.109469  0.091490  0.112792  0.097651  0.134221  0.163524
    ##           [,13]     [,14]     [,15]     [,16]     [,17]     [,18]
    ## C0    40.391425 40.102044 35.577207 50.117552 50.112882 47.827628
    ## C2    13.642838  7.941639  9.553459  6.372541  9.928400  8.704867
    ## C3     0.473142  0.475303  0.529805  0.586893  0.726657  0.611884
    ## C4     0.892873  0.821234  0.761605  1.363891  1.699377  0.984800
    ## C14:1  0.192286  0.156558  0.222255  0.147061  0.163120  0.315577
    ## C16    0.171114  0.141277  0.201599  0.091688  0.109650  0.130592
    ##           [,19]     [,20]     [,21]     [,22]     [,23]     [,24]  [,25]
    ## C0    31.749719 39.545500 38.466395 32.754739 46.383370 34.533526 45.400
    ## C2     4.524604  4.751794  5.632226  4.918422  5.630013  4.600135  4.760
    ## C3     0.417636  0.392226  0.301366  0.412230  0.564245  0.437928  0.571
    ## C4     0.522774  0.453093  0.171522  0.141726  0.190074  0.169483  0.218
    ## C14:1  0.137385  0.116437  0.080100  0.113568  0.106486  0.131095  0.085
    ## C16    0.069759  0.084442  0.093390  0.101643  0.082133  0.108082  0.122
    ##        [,26]  [,27]  [,28]  [,29]  [,30]  [,31]  [,32]  [,33]  [,34]
    ## C0    56.600 70.300 61.400 41.500 31.600 36.400 35.100 29.100 38.600
    ## C2     6.040  7.550  7.020  5.250  5.710  7.100  5.950  4.640  4.560
    ## C3     0.757  0.888  0.715  0.656  0.509  0.611  0.530  0.379  0.522
    ## C4     0.285  0.442  0.399  0.377  0.320  0.346  0.408  0.233  0.305
    ## C14:1  0.091  0.075  0.085  0.061  0.074  0.120  0.093  0.145  0.136
    ## C16    0.127  0.123  0.119  0.081  0.099  0.129  0.098  0.106  0.144
    ##        [,35]  [,36]  [,37]  [,38]  [,39]  [,40]  [,41]  [,42]  [,43]
    ## C0    32.500 32.000 35.900 27.600 54.300 52.200 60.000 17.900 37.800
    ## C2     3.510  5.730  7.210  6.750  8.570  9.370  8.200  2.980  4.220
    ## C3     0.562  0.412  0.659  0.462  0.826  0.805  0.646  0.122  0.547
    ## C4     0.412  0.338  0.441  0.392  0.600  0.675  0.725  0.128  0.358
    ## C14:1  0.091  0.135  0.140  0.150  0.124  0.130  0.106  0.050  0.059
    ## C16    0.138  0.124  0.145  0.104  0.150  0.144  0.132  0.062  0.071
    ##        [,44]  [,45]  [,46]  [,47]  [,48]  [,49]  [,50]  [,51]  [,52]
    ## C0    26.400 60.700 74.200 55.200 48.000 39.600 31.300 35.200 41.300
    ## C2     8.280 11.100 15.200 11.300  7.970  6.130  7.330  5.310  7.830
    ## C3     0.386  1.020  1.350  1.190  0.444  0.370  0.295  0.579  0.573
    ## C4     0.396  0.513  0.900  0.606  0.285  0.259  0.155  0.373  0.346
    ## C14:1  0.100  0.093  0.100  0.122  0.116  0.139  0.129  0.112  0.184
    ## C16    0.122  0.152  0.135  0.197  0.129  0.143  0.142  0.098  0.152
    ##        [,53]       [,54]       [,55]       [,56]     [,57]       [,58]
    ## C0    14.700 21.96500000 32.50800000 24.35900000 38.583000 27.97800000
    ## C2     3.540  4.16400000  6.96300000  3.10100000 17.261000  4.78900000
    ## C3     0.150  0.44878200  0.35359600  0.21777100  0.571344  0.37454400
    ## C4     0.264  0.49235900  1.05200000  0.20565800  0.514015  0.39871500
    ## C14:1  0.072  0.08649593  0.10909700  0.07763973  0.176367  0.09788839
    ## C16    0.075  0.04509888  0.09246439  0.08573237  0.173704  0.09900774
    ##             [,59]     [,60]       [,61]     [,62]
    ## C0    46.50700000 42.375000 34.42700000 37.334000
    ## C2     4.58300000  6.890000  7.87800000 11.554000
    ## C3     0.56977700  0.462423  0.34376400  0.430095
    ## C4     0.56726300  0.380776  0.29406000  0.453931
    ## C14:1  0.07387225  0.126970  0.09835858  0.154291
    ## C16    0.08533394  0.203613  0.14406100  0.166736

Access metabolite annotation

``` r
rowData(CPI_metabolomics)
```

    ## DataFrame with 134 rows and 8 columns
    ##              id          class       assay    length n_c22_isobars
    ##     <character>    <character> <character> <integer>     <integer>
    ## 1            C0 acylcarnitines        p180        NA            NA
    ## 2            C2 acylcarnitines        p180        NA            NA
    ## 3            C3 acylcarnitines        p180        NA            NA
    ## 4            C4 acylcarnitines        p180        NA            NA
    ## 5         C14:1 acylcarnitines        p180        NA            NA
    ## ...         ...            ...         ...       ...           ...
    ## 130    SM C24:0  sphingolipids        p180        24             2
    ## 131    SM C24:1  sphingolipids        p180        24             2
    ## 132    SM C26:0  sphingolipids        p180        26             3
    ## 133    SM C26:1  sphingolipids        p180        26             1
    ## 134          H1        hexoses        p180        NA            NA
    ##     n_isobars  fraction   bond_type
    ##     <integer> <numeric> <character>
    ## 1          NA        NA          NA
    ## 2          NA        NA          NA
    ## 3          NA        NA          NA
    ## 4          NA        NA          NA
    ## 5          NA        NA          NA
    ## ...       ...       ...         ...
    ## 130         2         1          NA
    ## 131         2         1          NA
    ## 132         3         1          NA
    ## 133         1         1          NA
    ## 134        NA        NA          NA

Access clinical annotation

``` r
colData(CPI_metabolomics)
```

    ## DataFrame with 62 rows and 3 columns
    ##               X      Gender       Age
    ##     <character> <character> <integer>
    ## 1       sample1           m        66
    ## 2       sample2           m        66
    ## 3       sample3           m        55
    ## 4       sample4           m        55
    ## 5       sample5           m        55
    ## ...         ...         ...       ...
    ## 58     sample58           m        58
    ## 59     sample59           m        55
    ## 60     sample60           m        65
    ## 61     sample61           m        64
    ## 62     sample62           m        77

Session information
===================

``` r
sessionInfo()
```

    ## R version 3.4.2 (2017-09-28)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS  10.14.2
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ## [1] SummarizedExperiment_1.8.1 DelayedArray_0.4.1        
    ## [3] matrixStats_0.54.0         Biobase_2.38.0            
    ## [5] GenomicRanges_1.30.1       GenomeInfoDb_1.14.0       
    ## [7] IRanges_2.12.0             S4Vectors_0.16.0          
    ## [9] BiocGenerics_0.24.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.19           knitr_1.20             XVector_0.18.0        
    ##  [4] magrittr_1.5           zlibbioc_1.24.0        lattice_0.20-35       
    ##  [7] stringr_1.3.1          tools_3.4.2            grid_3.4.2            
    ## [10] htmltools_0.3.6        yaml_2.2.0             rprojroot_1.3-2       
    ## [13] digest_0.6.18          Matrix_1.2-12          GenomeInfoDbData_1.0.0
    ## [16] bitops_1.0-6           RCurl_1.95-4.10        evaluate_0.11         
    ## [19] rmarkdown_1.10         stringi_1.2.4          compiler_3.4.2        
    ## [22] backports_1.1.2
