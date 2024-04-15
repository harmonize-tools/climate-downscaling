# climate variable (the same code is need as variable name of the netcdf file and in the file name)
  var_name = 't2m' # 2m temperature
 
# reference period, forecast issue date and leadtimes
  reference_period <- c(1996:2015)
  forecast_issue_date <- '2024-04'
  leadtimes <- indices(1:3)

# configuration of sdate_hcst (array containing the initialisation dates of the reference period)
  sdate_hcst <- paste0(reference_period, substr(forecast_issue_date,6,7))

# configuration of sdate_fcst (array containing the initialisation dates of the forecast)
  sdate_hcst <- paste0(substr(forecast_issue_date,1,4), substr(forecast_issue_date,6,7))
                           
# SELECTION of the region boundaries; the sample data is prepared for reginos inside lons (-85, -25) and lats (-15, 25)     
  lons.min <- -40       
  lons.max <- -32     
  lats.min <- -11   
  lats.max <- -3   

# path where to find the sample data
  exp_path <- paste0('./sample_data/ecmwf51/$var$_$sdate$01.nc')  
  obs_path <- paste0('./sample_data/era5land/$var$_$date$.nc')  
  
# Load and prepare hindcast data (seasonal prediction in the reference period)
  exp <- startR::Start(
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

# extract dates and coordinates

# Load Obs
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
# visualize Climatology (3-month avg, ensemble mean) of exp data and obs data and their resolutions

# Load forecast (3 months)
# visualize 

# downscale hindcast and assess quality
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

# Visualize metric quality assessment
# select final metric

# downscale forecast (3 months)
# Visualize raw forecast vs calibrated downscaled forecast
