### load libraries

#### dummy example
    exp <- rnorm(500)   
    dim(exp) <- c(member = 5, lat = 4, lon = 5, sdate = 5) 
    exp_lons <- 1:5 
    exp_lats <- 1:4 
    obs <- rnorm(900) 
    dim(obs) <- c(lat = 12, lon = 15, sdate = 5) 
    obs_lons <- seq(1,5, 4/14) 
    obs_lats <- seq(1,4, 3/11) 
    res_grid <- Intbc(exp = exp, obs = obs, exp_lats = exp_lats, exp_lons = exp_lons, obs_lats = obs_lats, 
              obs_lons = obs_lons, target_grid = 'r1280x640', bc_method = 'evmos', int_method = 'conservative')


    obs.point <- rnorm(20)
    dim(obs.point) <- c(lat = 2, lon = 2, sdate = 5)
    obs.point_lons <- c(1.4, 2)
    obs.point_lats <- c(3.2, 3)
    res_point <- Intbc(exp = exp, obs = obs.point, exp_lats = exp_lats, exp_lons = exp_lons, obs_lats = obs_lats, obs_lons = obs_lons,
                        target_grid = 'r1280x640', bc_method = 'evmos', int_method = 'dis',
                        points = list(lat = obs.point_lats, lon = obs.point_lons), method_point_interp = "nearest",
                        source_file_exp = './sample_data/ecmwf51/t2m_19960401.nc',
                        source_file_obs = './sample_data/era5land/t2m_199604.nc')
