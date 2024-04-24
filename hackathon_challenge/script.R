# reference period, forecast issue date and leadtimes
reference_period <- c(2002:2015)
forecast_issue_date <- '2024-04'
leadtimes <- 1

# configuration of sdate_hcst (array containing the initialisation dates of the reference period)
sdate_hcst <- paste0(reference_period, substr(forecast_issue_date,6,7))

# configuration of sdate_fcst (array containing the initialisation dates of the forecast)
sdate_fcst <- paste0(substr(forecast_issue_date,1,4), substr(forecast_issue_date,6,7))

# path where to find the sample data
exp_path <- paste0('./sample_data/ecmwf51/$var$_$sdate$01.nc')  
obs_path <- paste0('./sample_data/era5land/$var$_$date$.nc')

exp <- startR::Start(
  dat = exp_path,
  var = var_name,
  sdate = sdate_hcst,
  ensemble = 'all',
  time = leadtimes,
  latitude = values(list(latmin, latmax)),
  latitude_reorder = Sort(decreasing = T),
  longitude = values(list(lonmin, lonmax)),
  longitude_reorder = CircularSort(-180,180),
  synonims = list(latitude = c('lat', 'latitude'),
                  longitude = c('lon', 'longitude'),
                  ensemble=c('member','ensemble','number')),
  return_vars = list(latitude = 'dat',
                     longitude = 'dat',
                     time = 'sdate'),
  retrieve = TRUE)

points <- list(latitude = c(-3, -1), longitude = c(-73.3,-71))
obs_point <- point_data
down_points<- Intbc(exp = exp, obs = obs_point, exp_lats = attr(exp, "Variables")$dat1$lat,
                    exp_lons = attr(exp, "Variables")$dat1$lon, obs_lats = points$lat,
                    obs_lons = points$lon, points = points, method_point_interp = 'nearest',
                    source_file_exp = attr(exp, "Files")[1,1,1], source_file_obs = NULL,
                    target_grid = 'r4000x4000', 
                    bc_method = 'bias', ncores = 4, lon_dim = 'longitude', lat_dim = 'latitude', member_dim = 'ensemble')
