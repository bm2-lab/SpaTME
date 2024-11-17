# SpaLinker

## Overview

![Overview](inst/Figure_overview.jpg)
SpaLinker is an integrative analytical framework designed for TME deciphering and phenotype linking using spot-resolution ST data. It consists of three functional parts (Overview): (1) Characterization of the TME with molecular and cellular features; (2) Recognition of spatial architectures; and (3) Linking with clinical phenotypes.

## Installation

SpaLinker has been built and tested with R >= 4.3.2. Specific package dependencies are defined in the package DESCRIPTION.
```
# Install devtools, if necessary
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("bm2-lab/SpaLinker")
```

## Demonstration

For examples of typical SpaLinker usage, please see our [package vignette](https://bm2-lab.github.io/SpaLinker/vignettes/SpaLinker_tutorial) and [download](https://zenodo.org/uploads/13120480) the test data.

