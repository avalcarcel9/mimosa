FROM mimosa_base

RUN Rscript -e "chooseCRANmirror(graphics=FALSE, ind=56); \
                install.packages('https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz', repos=NULL, type='source'); \
                install.packages('https://cran.r-project.org/src/contrib/Archive/rlist/rlist_0.4.6.tar.gz', repos=NULL, type='source'); \
                source('https://neuroconductor.org/neurocLite.R'); neuro_install(c('ANTsRCore', 'extrantsr', 'mimosa'))"

