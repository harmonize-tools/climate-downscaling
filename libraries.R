# Binder:
  install.packages("./extra_libs/ncdf4_1.22.tar.gz",repos=NULL, type="source")
  pacman::p_load(ncdf4, startR, s2dv, CSTools, easyVerification, multiApply, ClimProjDiags, plyr, nnet, FNN, ecmwfr, devtools, lubridate)

# Local: 
  if(!require('pacman')){install.packages('pacman')}
  pacman::p_load(ncdf4, startR, s2dv, CSTools, easyVerification, multiApply, ClimProjDiags, plyr, nnet, FNN, ecmwfr, devtools, lubridate)

 ## alternatively:
    library(ncdf4)
    library(startR)
    library(s2dv)
    library(CSTools)
    library(easyVerification)
    library(multiApply)
    library(ClimProjDiags)
    library(plyr)
    library(nnet)
    library(FNN)
    library(ecmwfr)
    library(devtools)
    library(lubridate)

# source public code CSDownscale package (run either in Binder or Local):
  source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Analogs.R")
  source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Interpolation.R")
  source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Intbc.R")
  source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Intlr.R")
  source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/LogisticReg.R")
  source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Utils.R")

