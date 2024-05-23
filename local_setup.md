## Set up to run the tutorial locally
It is also possible and convenient to run the tutorial locally to reduce the potential issues due to connectivity. However, the installation of all the required packages is complex unless your system is Linux or you have access to a Linux VirtualBox. If you use Linux in your local computer, we recommend that you follow all this steps to ensure that you will be able to follow the tutorial locally. For other operating systems you can check the following links with information on how to install CDO in [Windows Systems](https://code.mpimet.mpg.de/projects/cdo/wiki/Win32) and [MacOS Platform](https://code.mpimet.mpg.de/projects/cdo/wiki/MacOS_Platform).

1. Check if CDO (Climate Data Operators) is installed: in the terminal type cdo -V to check the CDO version; if the command is not recognised, install CDO with sudo apt-get install cdo, then check again cdo -V (you should get some information about the version and system)

2. Check if you have a netcdf library configuration: in the terminal type which nc-config; if nothing appears, you need to install the configuration with sudo apt-get install libnetcdf-dev; check which nc-config again and a path information should appear (e.g. /usr/bin/nc-config or /usr/local/bin/nc-config)

3. If you don’t have it already, install R and the RStudio: https://posit.co/download/rstudio-desktop/

4. Open an R session or RStudio and install the necessary packages with the following lines:
      ```
      if(!require('pacman')){install.packages('pacman')}
      pacman::p_load(ncdf4, startR, s2dv, CSTools, easyVerification, multiApply, ClimProjDiags, plyr, nnet, FNN, ecmwfr, devtools, lubridate)
      ```

5. Get the sample data for the tutorial:
      - Create a directory where you would like to work
      - Click here to download the sample data for the tutorial (it will be automatically saved in your downloads folder)
      - Uncompress the downloaded file and save it in the directory where you want to work with the name sample_data and keeping the same structure (sample_data should contain ecmwf51 with 21 files and era5land with 60 files)
