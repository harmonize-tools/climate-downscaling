# Downscaling of Climate Data
Tutorial to spatially downscale climate data in the hotspots of the [HARMONIZE project](https://www.harmonize-tools.org/)

<!-- <img src="images/harmonize_logo.png" height="70"/> <img src="images/bsc_logo.png" height="70"/> -->

Click the following icon to run the tutorial with binder *(note that the binder session can take some time to load, about 15 to 30min, thank you for your patience)*: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/harmonize-tools/climate-downscaling/HEAD) 

<img src="images/csdownscale_logo.png" width="200"/> [CSDownscale](https://earth.bsc.es/gitlab/es/csdownscale)

Downscaling allows for transferring climate information from coarse to fine grids. In this way, seasonal predictions, which are usually delivered in coarse grids, can be output refined to improve their value. It is worth noting that downscaling does not necessarily increase the overall skill (quality) of seasonal forecasts. Instead, it provides a more detailed spatial field of climate variables like temperature, precipitation or surface winds. A wide variety of downscaling methods do exist. Some statistical methods have been coded and included in the CSDownscale R package.

- **Interpolation**: Included in Interpolation(). Regrid of a coarse-scale grid into a fine-scale grid, or interpolate model data into a point location. Different interpolation methods, based on different mathematical approaches, can be applied: conservative, nearest neighbour, bilinear or bicubic. Does not rely on any data for training.
- **Interpolation plus bias adjustment**: Included in Intbc(). interpolate model data into a fine-scale grid or point location. Later, a bias adjustment of the interpolated values is performed. Bias adjustment techniques include simple bias correction, calibration or quantile mapping.
- **Interpolation plus linear regression**: Included in Intlr(..., method = 'basic'). Firstly, model data is interpolated into a fine-scale grid or point location. Later, a linear-regression with the interpolated values is fitted using high-res observations as predictands, and then applied with model data to correct the interpolated values.
- **Interpolation with large-scale predictors**: Included in Intlr(..., method = 'large-scale'). Firstly, model data is interpolated into a fine-scale grid or point location. Later, a linear-regression with large-scale predictors from the same model (e.g. teleconnection indices) is fitted using high-res observations as predictands. Finally, the linear-regression is applied with model data to correct the interpolated values.
- **Stencil**: Included in Intlr(..., method = '9nn'). A linear-regression with the nine nearest neighbours is fitted using high-res observations as predictands. Instead of constructing a regression model using all the nine predictors, principal component analysis is applied to the data of neighbouring grids to reduce the dimension of the predictors. The linear regression model is then built using the principal components that explain 95% of the variance. The '9nn' method does not require a pre-interpolation process.
- **Analogs**: The analogs function determines the N best analogs based on Euclidian distance, distance correlation, or Spearman's correlation metrics. To downscale a local-scale variable, either the variable itself (i.e., Model Output Statistics (MOS) approach) or another large-scale variable (i.e., Perfect Prognosis (PP) approach) can be utilized as the predictor. In the first scenario, analogs are examined between the observation and model data of the same local-scale variable. In the latter scenario, the function identifies the day in the observation data that closely resembles the large-scale pattern of interest in the model. When it identifies the date of the best analog, the function extracts the corresponding local-scale variable for that day from the observation of the local scale variable. The used local-scale and large-scale variables can be retrieved from independent regions.
- **Logistic regression**: Included in LogisticReg(). this method uses a sigmoid function to relate ensemble mean anomalies of the large-scale forecasts directly to probabilities of observing above normal/normal/below normal conditions at the local scale. Therefore, it does not produce an ensemble of forecasts but rather their associated probabilities. It is a statistical method with few parameters to train, and only benefits from local information, but it has shown good performance.
