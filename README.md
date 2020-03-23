# tfa
tfa is an R package for performing factor analyses of temporal or ancient DNA samples, correcting individual scores for the effect of allele frequency drift through time. The package implements a fast algorithm that computes $K$ corrected factors, and provides estimates of ancestry proportions for a target individual or a target population given a set of source populations. See a short [Overview](https://bcm-uga.github.io/tfa/articles/main-vignette.html).

## Installation

Install the latest version from github (requires [devtools](https://github.com/hadley/devtools)):
```R
# install.packages("devtools")
devtools::install_github("bcm-uga/tfa")
```

## References

- Olivier François, Séverine Liégeois, Benjamin Demaille, Flora Jay. (2019). Inference of population genetic structure from temporal samples of DNA. bioRxiv 801324. [doi:10.1101/801324](https://doi.org/10.1101/801324).

- Olivier François, Flora Jay. (2020). Factor analysis of ancient population genomic samples. Under review.
