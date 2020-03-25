#!/bin/bash

docker run -it -v $(pwd):/src --rm pennsive/mimosa:latest bash -c "Rscript -e "'"'"install.packages('pkgdown'); devtools::install_github('r-lib/pkgdown'); pkgdown::build_site()"'"'""