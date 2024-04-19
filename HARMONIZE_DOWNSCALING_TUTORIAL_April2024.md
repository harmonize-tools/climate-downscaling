# Step 1: load the necessary libraries 
Check [libraries.R](https://github.com/harmonize-tools/climate-downscaling/blob/main/libraries.R) for the specific code and/or [SetUp.md](https://github.com/harmonize-tools/climate-downscaling/blob/main/SetUp.md) for more information.

# Step 2: define parameters
Run the following lines of code to create the necessary parameters:

```
#climate variable (the same code is need as variable name of the netcdf file and in the file name)
 var_name = 't2m' # 2m temperature
 
# reference period, forecast issue date and leadtimes
  reference_period <- c(1996:2015)
  forecast_issue_date <- '2024-04'
  leadtimes <- indices(1:3)

# configuration of sdate_hcst (array containing the initialisation dates of the reference period)
  sdate_hcst <- paste0(reference_period, substr(forecast_issue_date,6,7))

# configuration of sdate_fcst (array containing the initialisation dates of the forecast)
  sdate_hcst <- paste0(substr(forecast_issue_date,1,4), substr(forecast_issue_date,6,7))

# path where to find the sample data
  exp_path <- paste0('./sample_data/ecmwf51/$var$_$sdate$01.nc')  
  obs_path <- paste0('./sample_data/era5land/$var$_$date$.nc')
```
                           
# Step 3: SELECTION of the region boundaries
The sample data is prepared for reginos inside lons (-85, -25) and lats (-15, 25). See [Hotspots.md](https://github.com/harmonize-tools/climate-downscaling/blob/main/Hotspots.md) for more information.

Run the following lines of code after modifying them if you want to select a different region:

```
  lons.min <- -78    
  lons.max <- -72  
  lats.min <- 1.5   
  lats.max <- 7.5  
```
  
# Step 4: Load data into the session
## Load and prepare hindcast data (seasonal prediction in the reference period)
```
  hcst <- startR::Start(
          dat = exp_path,
          var = var_name,
          sdate = sdate_hcst,
          ensemble = 'all',
          time = leadtimes,
          latitude = values(list(lats.min, lats.max)),
          latitude_reorder = Sort(decreasing = T),
          longitude = values(list(lons.min, lons.max)),
          longitude_reorder = CircularSort(-180,180),
          synonims = list(latitude = c('lat', 'latitude'),
                          longitude = c('lon', 'longitude'),
                          ensemble=c('member','ensemble','number')),
          return_vars = list(latitude = 'dat',
                              longitude = 'dat',
                              time = 'sdate'),
          retrieve = TRUE)

  # transform the units (from Kelvin to Celsius)
  if (attr(hcst, "Variables")$common[[2]]$units == 'K'){
      hcst <- hcst - 273.15
      attr(hcst, "Variables")$common[[2]]$units <- 'C'
  }
  # extract dates and coordinates
  dates_hcst <- attr(hcst, 'Variables')$common$time
  lats_hcst <- attr(hcst, "Variables")$dat1$latitude
  lons_hcst <- attr(hcst, "Variables")$dat1$longitude
```

## Load and prepare reanalysis data to be used as reference
```
  # Here, the dates of the previously loaded hindcast will be used
  # to retrieve the reanalysis. They will be the same as the hindcast times.
  dates_file <- format(dates_hcst, '%Y%m') #Giving dates format
  dim(dates_file) <- c(sdate = length(sdate_hcst), time = length(leadtimes))#Specifying the dimensions 
  
  obs <- Start(# load observational (reanalysis) data
      dat = obs_path,
      var = var_name,
      date = dates_file,
      latitude = values(list(lats.min,lats.max)),
      latitude_reorder = Sort(decreasing = T),
      longitude = values(list(lons.min, lons.max)),
      longitude_reorder = CircularSort(-180,180),
      synonims = list(longitude = c('lon', 'longitude'),
                      latitude = c('lat', 'latitude')),
      split_multiselected_dims = TRUE,
      return_vars = list(time = 'date',
                          latitude = 'dat',
                          longitude = 'dat'),
      retrieve = TRUE)

  # transform the units (from Kelvin to Celsius)
  if (attr(obs, "Variables")$common[[2]]$units == 'K'){
      obs <- obs - 273.15
      attr(obs, "Variables")$common[[2]]$units <- 'C'
  }
  # extract dates and coordinates
  dates_obs <- attr(obs, 'Variables')$common$time
  lats_obs <- attr(obs, "Variables")$dat1$latitude
  lons_obs <- attr(obs, "Variables")$dat1$longitude 
```

## Load the forecast (a seasonal prediction into the future)
```
  fcst <- startR::Start(
      dat = exp_path,
      var = var_name,
      sdate = sdate_fcst,
      ensemble = 'all',
      time = leadtimes,
      latitude = values(list(lats.min, lats.max)),
      latitude_reorder = Sort(decreasing = T),
      longitude = values(list(lons.min, lons.max)),
      longitude_reorder = CircularSort(-180,180),
      synonims = list(latitude = c('lat', 'latitude'),
                      longitude = c('lon', 'longitude'),
                      ensemble=c('member','ensemble','number')),
      return_vars = list(latitude = 'dat',
                          longitude = 'dat',
                          time = 'sdate'),
      retrieve = TRUE)

  # transform the units (from Kelvin to Celsius)
  if (attr(fcst, "Variables")$common[[2]]$units == 'K'){
      fcst <- fcst - 273.15
      attr(fcst, "Variables")$common[[2]]$units <- 'C'
  }
  # extract dates and coordinates
  dates_fcst <- attr(fcst, 'Variables')$common$time
  lats_fcst <- attr(fcst, "Variables")$dat1$latitude
  lons_fcst <- attr(fcst, "Variables")$dat1$longitude 
```

## Calculate ensemble mean, climatologies, anomalies and seasonal averages for visualisations
```  
  hcst.ensemble_mean <- MeanDims(hcst, dim = 'ensemble', na.rm = TRUE)
  fcst.ensemble_mean <- MeanDims(fcst, dim = 'ensemble', na.rm = TRUE)

  hcst.clim <- MeanDims(hcst.ensemble_mean, dim = 'sdate', na.rm = TRUE)
  obs.clim <- MeanDims(obs, dim = 'sdate', na.rm = TRUE)

  hcst.anom <- Ano(data = hcst, clim = hcst.clim)
  fcst.anom <- Ano(data = fcst, clim = hcst.clim)
  obs.anom <- Ano(data = obs, clim = obs.clim)
  
  hcst.clim_season <- MeanDims(hcst.clim, dim = 'time', na.rm = TRUE)
  obs.clim_season <- MeanDims(obs.clim, dim = 'time', na.rm = TRUE)
  fcst.season_av <- MeanDims(fcst, dim = 'time', na.rm = TRUE)
```

# Step 3: Create visualise some plots of the loaded raw data

# Visualize metric quality assessment
# select final metric

# downscale forecast (3 months)
# Visualize raw forecast vs calibrated downscaled forecast
