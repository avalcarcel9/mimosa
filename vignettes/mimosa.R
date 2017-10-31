## ---- eval = FALSE-------------------------------------------------------
#  source("https://neuroconductor.org/neurocLite.R")
#  neuro_install("mimosa")

## ---- eval = FALSE-------------------------------------------------------
#  devtools::install_github('avalcarcel9/mimosa')

## ---- warning = FALSE, message = FALSE-----------------------------------
library(neurobase)
library(mimosa)
# library(data.table)
library(dplyr)
library(oasis)
library(fslr)
# library(scales)
# library(tidyverse)

## ----setup, eval = FALSE, echo = FALSE-----------------------------------
#  knitr::opts_knit$set(eval = FALSE)

## ---- echo = FALSE-------------------------------------------------------
info = Sys.info()
user = info[["user"]]
if (grepl("muschel", user )) {
  data_dir = "~/Desktop/"
} else {
  data_dir = "/Users/alval/Documents"
}

## ------------------------------------------------------------------------
# Note these paths will be to where you 
# have stored the data or to your own data
train_dir = file.path(data_dir, "training01")

T1_files = list.files(path = 
                        file.path(train_dir,
                                  "preprocessed"), 
                      pattern = "mprage_pp[.]nii", 
                      full.names = TRUE)
T2_files = list.files(path = 
                        file.path(train_dir,
                                  "preprocessed"), 
                      pattern = "t2_pp[.]nii", 
                      full.names = TRUE)
FLAIR_files = list.files(path = 
                           file.path(train_dir,
                                     "preprocessed"), 
                         pattern = "flair_pp[.]nii", 
                         full.names = TRUE)
PD_files = list.files(path = 
                        file.path(train_dir,
                                  "preprocessed"), 
                      pattern = "pd_pp[.]nii", 
                      full.names = TRUE)
GS_files = list.files(path = 
                        file.path(train_dir,
                                  "masks"), 
                      pattern = "mask2[.]nii", 
                      full.names = TRUE)
filepaths = data.frame(T1 = T1_files, T2 = T2_files, FLAIR = FLAIR_files, PD = PD_files, GS = GS_files, stringsAsFactors = FALSE)
ss = strsplit(nii.stub(filepaths$T1), split = "_")
filepaths$visit_id = sapply(ss, function(x) x[2])
filepaths

