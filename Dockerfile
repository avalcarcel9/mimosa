FROM pennsive/r-env:base

# install anything not already in pennsive/r-env:base
RUN r -e "install.packages('rlist')"
RUN r -e "devtools::install_github('avalcarcel9/mimosa', dependencies = FALSE)"

WORKDIR /src
COPY . .
ENTRYPOINT []
CMD bash
