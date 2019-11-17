HawkesFit tutorial: simulation, fitting and extension
================

## Introduction

This package is designed for simulating and fitting the Hawkes processes
and the HawkesN processes with several options of kernel functions.
Currently, it assumes univariate processes without background event
rates. Prior knowledge about the models is assumed in the following
tutorial and please refer to \[1\] and \[2\] for details about the
models.

## Simulation

``` r
lapply(seq(2), function(i) generate_Hawkes_event_series(c(K = 0.8, theta = 1), model_type = 'EXP'))
```

    ## [[1]]
    ##      magnitude      time
    ## 1 10000.000000 0.0000000
    ## 2     1.035351 0.2532774
    ## 
    ## [[2]]
    ##   magnitude time
    ## 1     10000    0

## Including Code

You can include R code in the document as follows:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Reference

\[1\] Rizoiu, Marian-Andrei, et al. “Hawkes processes for events in
social media.” Frontiers of Multimedia Research. Association for
Computing Machinery and Morgan & Claypool, 2017.

\[2\] Rizoiu, Marian-Andrei, et al. “SIR-Hawkes: Linking epidemic models
and Hawkes processes to model diffusions in finite populations.”
Proceedings of the 2018 World Wide Web Conference. International World
Wide Web Conferences Steering Committee, 2018. You can also embed plots,
for example:
