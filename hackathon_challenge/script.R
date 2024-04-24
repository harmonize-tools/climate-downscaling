
point_data_1 <- as.array(read.csv("./quistocoha.csv") [,3])
point_data_2 <- as.array(read.csv("./el_paujil.csv") [,3])
dim(point_data_1)  <- dim(point_data_2 )<- c(length(point_data_1))
point_data <- abind (point_data_1, point_data_2, along=2)
names(dim(point_data)) <- c("sdate", "location")
dim(point_data) <- c(dim(point_data), dat = 1, time = 1, var = 1 ) 
#target_grid <- '/esarchive/recon/ecmwf/era5/monthly_mean/tas_f1h/tas_200002.nc'
lonmin <- -90
lonmax <- -30
latmin <- -15
latmax <- 0
#sdates <- c('20000201','20010201','20020201','20030201','20040201','20050201','20060201','20070201')
sdates <- format(ymd("20000201") + rep(years(0:13), each=1),"%Y%m%d")
exp <- startR::Start(dat = '/esarchive/exp/ecmwf/system5c3s/monthly_mean/$var$_f6h/$var$_$sdate$.nc',
                     var = 'tas', time = indices(1), member = 'all', sdate = sdates,
                     lat = values(list(latmin, latmax)), lat_reorder = Sort(decreasing = FALSE),
                     lon = values(list(lonmin, lonmax)), lon_reorder = CircularSort(0, 360),
                     synonims = list(var = c('var','variable'), lon = c('lon', 'longitude'),
                     lat = c('lat', 'latitude'), member = c('member','ensemble')),
                     return_vars = list(lat = 'dat', lon = 'dat'),
                     num_procs = 1, retrieve = TRUE)
points <- list(lat = c(-3, -1), lon = c(360-73.3,360-71))
obs_point <- point_data
down_points<- Intbc(exp = exp, obs = obs_point, exp_lats = attr(exp, "Variables")$dat1$lat,
                          exp_lons = attr(exp, "Variables")$dat1$lon, obs_lats = points$lat,
                          obs_lons = points$lon, points = points, method_point_interp = 'nearest',
                          source_file_exp = attr(exp, "Files")[1,1,1], source_file_obs = NULL,
                          target_grid = target_grid, bc_method = 'bias', ncores = 4)
