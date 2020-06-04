# tfa: factor analysis of temporal samples

## Summary
The R package `tfa` implements a factor analysis algorithm for temporal DNA or ancient DNA samples, ajusting individual scores for the effect of allele frequency drift through time. The adjusted scores provide geometric representations of temporal samples consistent with estimates of ancestry proportions and have  interpretation similar to principal component analysis. Based on the adjusted factors, the program can also estimate ancestry proportions for a target population or a subset of target individuals given specified source populations. 


## Package overview

The main functions of the package are 
- ancestry_coefficients(): Compute ancestry coefficients from specified source populations
- coverage_adjust(): Adjust ancient genotypes for coverage bias
- tfa(): Factor analysis of temporal or ancient DNA samples

Check the following link to access the package website and a short [tutorial](https://bcm-uga.github.io/tfa/articles/tfa-vignette.html)

## Installation

Install the latest version from github (requires [devtools](https://github.com/hadley/devtools)):
```R
# install.packages("devtools")
devtools::install_github("bcm-uga/tfa")
```

## References

- Olivier François, Séverine Liégeois, Benjamin Demaille, Flora Jay. (2019). Inference of population genetic structure from temporal samples of DNA. bioRxiv 801324. [doi:10.1101/801324](https://doi.org/10.1101/801324).

- Olivier François, Flora Jay. (2020). Factor analysis of ancient population genomic samples. Under review.
