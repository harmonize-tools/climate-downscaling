if(!require('pacman')){install.packages('pacman')}
pacman::p_load(ncdf4, startR, s2dv, CSTools, easyVerification, multiApply, ClimProjDiags, plyr, nnet, FNN, ecmwfr, devtools, lubridate)

source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Analogs.R")
source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Interpolation.R")
source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Intbc.R")
source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Intlr.R")
source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/LogisticReg.R")
source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Utils.R")
