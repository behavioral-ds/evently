evently: simulation, fitting of Hawkes processes
![R-CMD-check](https://github.com/behavioral-ds/evently/workflows/R-CMD-check/badge.svg)
================

## Introduction

This package is designed for simulating and fitting the Hawkes processes
and the HawkesN processes with several options of kernel functions.
Currently, it assumes univariate processes without background event
rates. Prior knowledge about the models is assumed in the following
tutorial and please refer to `[1]` and `[2]` for details about the
models.

``` r
library(evently)
```

## Citation

Please consider citing the following paper if you find the package helpful to your research
```
@inproceedings{kong2021evently,
    address = {Jerusalem, Israel},
    author = {Kong, Quyu and Ram, Rohit and Rizoiu, Marian-Andrei},
    booktitle = {ACM International Conference on Web Search and Data Mining (WSDM), Demo},
    title = {{Evently: Modeling and Analyzing Reshare Cascades with Hawkes Processes}},
    year = {2021},
}
```


## Installation and dependencies

Several dependencies
([poweRlaw](https://cran.r-project.org/web/packages/poweRlaw/poweRlaw.pdf),
[AMPL](https://ampl.com/),
[Ipopt](https://www.coin-or.org/Ipopt/documentation/)) are required for
running this package. These dependencies will be installed automatically
by R or by following instructions upon package load.

Install the package by executing

``` r
if (!require('devtools')) install.packages('devtools')
devtools::install_github('behavioral-ds/evently')
```

## Simulating cascades

Let’s first simulate 100 event cascades of the **Hawkes process with an
exponential kernel function** (please refer to the [Available
models](#available-models) for models and their abbreviations in the
package) with a given parameter set, ![\\kappa = 0.9, \\theta
= 1](https://latex.codecogs.com/png.latex?%5Ckappa%20%3D%200.9%2C%20%5Ctheta%20%3D%201
"\\kappa = 0.9, \\theta = 1"). For each simulation, we only simulate
until 5 seconds. The resulted cascades are placed in a single `list`
where each cascade is a `data.frame`.

``` r
set.seed(4)
sim_no <- 100
data <- generate_series(par = c(K = 0.9, theta = 1), model_type = 'EXP', Tmax = 5, sim_no = sim_no)
# alternatively, `generate_series` also accepts a model class object
# e.g.
# model <- new_hawkes(par = c(K = 0.9, theta = 1), model_type = 'EXP')
# generate_series(model = model, Tmax = 5, sim_no = sim_no)

head(data[[1]])
```

    ##   magnitude      time
    ## 1         1 0.0000000
    ## 2         1 0.5941959
    ## 3         1 1.4712411
    ## 4         1 1.6105430
    ## 5         1 1.7855535
    ## 6         1 1.8883869

A simulated process is represented by a `data.frame` where each row is
an event. `time` indicates the event happening time, while `magnitude`
is the event mark information which is always 1 if `model_type` is an
unmarked model. In the context of retweet diffusion cascades, the first
row is the original tweet and all following events are its retweets.
`time` records the relative time (in second) of each retweet to the
original tweet and `magnitude` is the follows’ count of the user who
retweeted.

## Fitting a model on data

We can then fit on the cascades simulated in the previous section. After
providing the `data` and `model_type`, the fitting procedure will spawn
10 AMPL optimization procedures with different parameter
inistializations due to the non-convexity of some likelihood functions.
Among the 10 fitted model, the one giving the best likelihood value will
be returned. To make the fitting procedure faster, we can specify the
number of `cores` to be used for fitting them in
parallel.

``` r
fitted_model <- fit_series(data, model_type = 'EXP', observation_time = 5, cores = 10)

fitted_model
```

    ## - Model: EXP 
    ## - No. of cascades: 100 
    ## - init_par:
    ##   K 7.92e+00; theta 1.32e+00
    ## - par:
    ##   K 8.51e-01; theta 1.06e+00
    ## - Neg Log Likelihood: 285.488 
    ## - lower_bound:
    ##   K 1.00e-100; theta 1.00e-100
    ## - upper_bound:
    ##   K 1.00e+04; theta 3.00e+02
    ## - Convergence: 0

## Cascades from real data

We provide a function to help parse provided Tweet JSON objects to
cascades. Tweet JSON objects are crawled from the Twitter API and one
can refer to the
[documentation](https://developer.twitter.com/en/docs/tweets/data-dictionary/overview/tweet-object)
for more information.

In this package, we provide a collection of tweets about the same topic
for demonstration purpose. The file is in a
[`jsonl`](http://jsonlines.org/) format where each line of this file is
a Tweet JSON string. The content of the tweets have been encrypted due
to the Twitter API agreement. One can find the file at
`inst/extdata/tweets_anonymized.jsonl`. Let’s now extract cascades from
the
file.

``` r
filepath <- system.file('extdata', 'tweets_anonymized.jsonl', package = 'evently')
cascades <- parse_raw_tweets_to_cascades(filepath, progress = F)
print(cascades[seq(3)])
```

    ## [[1]]
    ##   magnitude time
    ## 1       337    0
    ## 
    ## [[2]]
    ##   magnitude time
    ## 2       130    0
    ## 
    ## [[3]]
    ##   magnitude time
    ## 3       251    0

As the file includes all raw tweets, it is natural to see that most
cascades extracted only have single events in it. We can further fit
Hawkes models on these cascades.

## Available models

There are 8 models available so far in this
package:

|                           Model                            | Abbreviation (model\_type) |                                                                                                                                                                  Intensity Function                                                                                                                                                                  |    Parameters    |
| :--------------------------------------------------------: | :------------------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :--------------: |
|     Hawkes process with an exponential kernel function     |            EXP             |                                                      ![\\kappa\\sum\_{t\_i \< t} \\theta e^{-\\theta (t-t\_i)}](https://latex.codecogs.com/png.latex?%5Ckappa%5Csum_%7Bt_i%20%3C%20t%7D%20%5Ctheta%20e%5E%7B-%5Ctheta%20%28t-t_i%29%7D "\\kappa\\sum_{t_i \< t} \\theta e^{-\\theta (t-t_i)}")                                                       |     K,theta      |
|      Hawkes process with a power-law kernel function       |             PL             |                                                      ![\\kappa\\sum\_{t\_i \< t} (t-t\_i + c)^{-(1+\\theta)}](https://latex.codecogs.com/png.latex?%5Ckappa%5Csum_%7Bt_i%20%3C%20t%7D%20%28t-t_i%20%2B%20c%29%5E%7B-%281%2B%5Ctheta%29%7D "\\kappa\\sum_{t_i \< t} (t-t_i + c)^{-(1+\\theta)}")                                                      |    K,c,theta     |
|    HawkesN process with an exponential kernel function     |            EXPN            |                         ![\\kappa\\frac{N-N\_t}{N}\\sum\_{t\_i \< t} \\theta e^{-\\theta (t-t\_i)}](https://latex.codecogs.com/png.latex?%5Ckappa%5Cfrac%7BN-N_t%7D%7BN%7D%5Csum_%7Bt_i%20%3C%20t%7D%20%5Ctheta%20e%5E%7B-%5Ctheta%20%28t-t_i%29%7D "\\kappa\\frac{N-N_t}{N}\\sum_{t_i \< t} \\theta e^{-\\theta (t-t_i)}")                          |    K,theta,N     |
|      HawkesN process with a power-law kernel function      |            PLN             |                         ![\\kappa\\frac{N-N\_t}{N}\\sum\_{t\_i \< t} (t-t\_i + c)^{-(1+\\theta)}](https://latex.codecogs.com/png.latex?%5Ckappa%5Cfrac%7BN-N_t%7D%7BN%7D%5Csum_%7Bt_i%20%3C%20t%7D%20%28t-t_i%20%2B%20c%29%5E%7B-%281%2B%5Ctheta%29%7D "\\kappa\\frac{N-N_t}{N}\\sum_{t_i \< t} (t-t_i + c)^{-(1+\\theta)}")                         |   K,c,theta,N    |
| Marked Hawkes process with an exponential kernel function  |            mEXP            |                              ![\\kappa\\sum\_{t\_i \< t} \\theta m\_i^{\\beta} e^{-\\theta (t-t\_i)}](https://latex.codecogs.com/png.latex?%5Ckappa%5Csum_%7Bt_i%20%3C%20t%7D%20%5Ctheta%20m_i%5E%7B%5Cbeta%7D%20e%5E%7B-%5Ctheta%20%28t-t_i%29%7D "\\kappa\\sum_{t_i \< t} \\theta m_i^{\\beta} e^{-\\theta (t-t_i)}")                              |   K,beta,theta   |
|   Marked Hawkes process with a power-law kernel function   |            mPL             |                             ![\\kappa\\sum\_{t\_i \< t} m\_i^{\\beta} (t-t\_i + c)^{-(1+\\theta)}](https://latex.codecogs.com/png.latex?%5Ckappa%5Csum_%7Bt_i%20%3C%20t%7D%20m_i%5E%7B%5Cbeta%7D%20%28t-t_i%20%2B%20c%29%5E%7B-%281%2B%5Ctheta%29%7D "\\kappa\\sum_{t_i \< t} m_i^{\\beta} (t-t_i + c)^{-(1+\\theta)}")                              |  K,beta,c,theta  |
| Marked HawkesN process with an exponential kernel function |           mEXPN            | ![\\kappa\\frac{N-N\_t}{N}\\sum\_{t\_i \< t} \\theta m\_i^{\\beta} e^{-\\theta (t-t\_i)}](https://latex.codecogs.com/png.latex?%5Ckappa%5Cfrac%7BN-N_t%7D%7BN%7D%5Csum_%7Bt_i%20%3C%20t%7D%20%5Ctheta%20m_i%5E%7B%5Cbeta%7D%20e%5E%7B-%5Ctheta%20%28t-t_i%29%7D "\\kappa\\frac{N-N_t}{N}\\sum_{t_i \< t} \\theta m_i^{\\beta} e^{-\\theta (t-t_i)}") |  K,beta,theta,N  |
|  Marked HawkesN process with a power-law kernel function   |            mPLN            |   ![\\kappa\\frac{N-N\_t}{N}\\sum\_{t\_i \< t} m\_i^{\\beta}(t-t\_i + c)^{-(1+\\theta)}](https://latex.codecogs.com/png.latex?%5Ckappa%5Cfrac%7BN-N_t%7D%7BN%7D%5Csum_%7Bt_i%20%3C%20t%7D%20m_i%5E%7B%5Cbeta%7D%28t-t_i%20%2B%20c%29%5E%7B-%281%2B%5Ctheta%29%7D "\\kappa\\frac{N-N_t}{N}\\sum_{t_i \< t} m_i^{\\beta}(t-t_i + c)^{-(1+\\theta)}")   | K,beta,c,theta,N |

## Acknowledgement

The development of this package is supported by the Green Policy grant
from the National Security College, Crawford School, ANU.

## License

Both dataset and code are distributed under the [Creative Commons
Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)
license](https://creativecommons.org/licenses/by-nc/4.0/). If you
require a different license, please contact us at <Quyu.Kong@anu.edu.au>
or <Marian-Andrei@rizoiu.eu>.

## Documentation

Please consult the package
[documentation](https://www.behavioral-ds.science/evently) for more details
and tutorials

## Reference

> `[1]` Rizoiu, M. A., Lee, Y., Mishra, S., & Xie, L. (2017, December).
> Hawkes processes for events in social media. In Frontiers of
> Multimedia Research (pp. 191-218). Association for Computing Machinery
> and Morgan & Claypool.  
> `[2]` Rizoiu, M. A., Mishra, S., Kong, Q., Carman, M., & Xie, L.
> (2018, April). SIR-Hawkes: Linking epidemic models and Hawkes
> processes to model diffusions in finite populations. In Proceedings of
> the 2018 World Wide Web Conference (pp. 419-428). International World
> Wide Web Conferences Steering Committee.  
> `[3]` Mishra, S., Rizoiu, M. A., & Xie, L. (2016, October). Feature
> driven and point process approaches for popularity prediction. In
> Proceedings of the 25th ACM International on Conference on Information
> and Knowledge Management (pp. 1069-1078). ACM.  
> `[4]` Kong, Q., Rizoiu, M. A., & Xie, L. (2019). Modeling Information
> Cascades with Self-exciting Processes via Generalized Epidemic Models.
> arXiv preprint arXiv:1910.05451.
