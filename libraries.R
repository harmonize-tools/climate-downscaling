# Binder:
  install.packages("./extra_libs/ncdf4_1.22.tar.gz",repos=NULL, type="source")
  library(startR, s2dv, CSTools, easyVerification, multiApply, ClimProjDiags, plyr, nnet, FNN, ncdf4, ecmwfr, devtools, lubridate)

# Local: 
  if(!require('pacman')){install.packages('pacman')}
  pacman::p_load(startR, s2dv, CSTools, easyVerification, multiApply, ClimProjDiags, plyr, nnet, FNN, ncdf4, ecmwfr, devtools, lubridate)

# source public code CSDownscale package (run either in Binder or Local):
  source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Analogs.R")
  source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Interpolation.R")
  source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Intbc.R")
  source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Intlr.R")
  source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/LogisticReg.R")
  source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Utils.R")
