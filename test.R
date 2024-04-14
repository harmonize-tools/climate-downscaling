#setwd("")

# variable selection 
  var_name = 't2m' 

# variable configuration (although this is repeated in download_seasonal_data.R)
  switch (tolower(var_name),
          "t2m"    = {var_longname <- '2m_temperature'},
          "tp"   = {var_longname <- 'total_precipitation'},
          "ssr"    = {var_longname <- 'surface_net_solar_radiation'},
          "msl"   = {var_longname <- 'mean_sea_level_pressure'})
  
# Selection of reference period, smonth (ini) and leadtimes (forecast_time)
  years <- 1996
  ini = '04'
  forecast_time <- indices(1:3)  

# configuration of sdate_forecast (array containing the initialisation dates)
  sdate_forecast <- paste0(years,ini)  
  
# Selection of the region boundaries
  lons.min <- -40       
  lons.max <- -32     
  lats.min <- -11   
  lats.max <- -3   

# path where the data will be stored SHOULD NOT BE CHANGED
  exp_path <- paste0('./sample_data/ecmwf51/$var$_$sdate$01.nc')  
  obs_path <- paste0('./sample_data/era5land/$var$_$date$.nc')  
  
# 4. Load and prepare seasonal data
  exp <- startR::Start(
    # Select the path of the forecast
    dat = exp_path,
    # Variable of interest
    var = var_name,
    # Start dates, years that we are using to calibrate the hindcast with reanalysis
    sdate = sdate_forecast,
    # Select all ensemble members
    ensemble = 'all',
    # Forecast time
    time = forecast_time,
    latitude = values(list(lats.min, lats.max)),
    latitude_reorder = Sort(decreasing = T),
    longitude = values(list(lons.min, lons.max)),
    # Reorder longitude points from [0,360] to [-180, 180]
    longitude_reorder = CircularSort(-180,180),
    synonims = list(latitude = c('lat', 'latitude'),
                    longitude = c('lon', 'longitude'),
                    ensemble=c('member','ensemble','number')),
    return_vars = list(latitude = 'dat',
                       longitude = 'dat',
                       time = 'sdate'),
    retrieve = TRUE)

  # Here, the dates that will be used to retrieve the reanalysis are obtained.
  # They will be the same as the hindcast times.
  dates <- attr(exp, 'Variables')$common$time
  
  #Giving dates format
  dates_file <- format(dates, '%Y%m') 
  
  #Specifying the dimensions
  dim(dates_file) <- c(sdate = length(sdate_forecast),
                       time = length(forecast_time))
  obs <- Start(# load observational (reanalysis) data
    dat = obs_path,
    # Variable of interest
    var = var_name,
    date = dates_file,
    # time = values(dates),
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
  exp_lats <- attr(exp, "Variables")$dat1$latitude
  exp_lons <- attr(exp, "Variables")$dat1$longitude
  obs_lats <- attr(obs, "Variables")$dat1$latitude
  obs_lons <- attr(obs, "Variables")$dat1$longitude
  clim_exp <- MeanDims(exp, dims = c('sdate', 'ensemble'), na.rm = TRUE)
  clim_obs <- MeanDims(obs, dims = c('sdate'), na.rm = TRUE)
  exp_anom <- Ano(data = exp, clim = clim_exp)
  obs_anom <- Ano(data = obs, clim = clim_obs)
  
  downscaled_hcst <- Intbc(exp = exp_anom, exp_lats = exp_lats, exp_lons = exp_lons,
                           obs = obs_anom, obs_lats = obs_lats, obs_lons = obs_lons,            
                           target_grid = 'r3600x1801', 
                           int_method = 'dis', 
                           lat_dim = 'latitude', 
                           member_dim = "ensemble",
                           bc_method = 'evmos', 
                           lon_dim = 'longitude', 
                           sdate_dim = 'sdate', 
                           ncores = 1)
