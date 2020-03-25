# tfa
The R package `tfa` implements a factor analysis algorithm for temporal DNA or ancient DNA (aDNA) samples, ajusting individual scores for the effect of allele frequency drift through time. Based on the adjusted factors, the program can estimate ancestry proportions for a target population or a subset of target individuals given specified source populations. Check the following link for a short [Overview](https://bcm-uga.github.io/tfa/articles/tfa-vignette.html) of the package.

## Installation

Install the latest version from github (requires [devtools](https://github.com/hadley/devtools)):
```R
# install.packages("devtools")
devtools::install_github("bcm-uga/tfa")
```

## References

- Olivier François, Séverine Liégeois, Benjamin Demaille, Flora Jay. (2019). Inference of population genetic structure from temporal samples of DNA. bioRxiv 801324. [doi:10.1101/801324](https://doi.org/10.1101/801324).

- Olivier François, Flora Jay. (2020). Factor analysis of ancient population genomic samples. Under review.
