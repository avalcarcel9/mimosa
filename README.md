# mimosa

mimosa is an R package based on the paper: MIMoSA: A Method for Inter-Modal Segmentation Analysis by Valcarcel et al. This package creates data structures necessary for training and testing and then allows the user to train a model and then apply the trained model to generate probability maps and predicted lesion segmentations.

## Installation

To install the package from neuroconductor, type:
```{r, eval = FALSE}
source("https://neuroconductor.org/neurocLite.R")
neuro_install("mimosa")
```

To get the latest development version from GitHub:

```{r, eval = FALSE}
devtools::install_github('avalcarcel9/mimosa')
```


[![Travis-CI Build Status](https://travis-ci.org/avalcarcel9/mimosa.svg?branch=master)](https://travis-ci.org/avalcarcel9/mimosa)

[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/avalcarcel9/mimosa?branch=master&svg=true)](https://ci.appveyor.com/project/avalcarcel9/mimosa)

[![Coverage Status](https://img.shields.io/coveralls/muschellij2/mimosa.svg)](https://coveralls.io/r/muschellij2/mimosa?branch=master)

## Vignette

For a full implementation of the methods please see our [vignette](https://github.com/avalcarcel9/mimosa/tree/master/vignettes).





