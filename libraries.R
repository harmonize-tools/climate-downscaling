# install required packages (if they are not yet installed) and load them:
install.and.load <- function(packagename){
  if (!require(packagename)){
    install.packages(packagename)
    library(packagename)
  } else {
    library(packagename)
  }
}
install.and.load('startR')
install.and.load('s2dv')
install.and.load('CSTools')
install.and.load('easyVerification')
install.and.load('multiApply')
install.and.load('ClimProjDiags')
install.and.load('plyr')
install.and.load('nnet')
install.and.load('FNN')
install.and.load('ncdf4')
install.and.load('ecmwfr')
install.and.load('devtools')
install.and.load('lubridate')

# source public code (not available as a package):
source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Analogs.R")
source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Interpolation.R")
source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Intbc.R")
source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Intlr.R")
source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/LogisticReg.R")
source("https://earth.bsc.es/gitlab/es/csdownscale/-/raw/master/R/Utils.R")

source('https://earth.bsc.es/gitlab/es/csdownscale/-/raw/dev-vignette/vignette/download_seasonal_data.R')

