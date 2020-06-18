## Overview

`rsplash` is the R implementation of the Simple process-led algorithms for simulating habitats (SPLASH v.2.0), which comprises robust formulations to compute energy and water fluxes. This R package, wrapping the C++ code, is intended to provide simulations either at site-scale or spatially-distributed, when the grid functionality is used. For reference, the code of the original v.1.0  (Davis et a., 2017) is hosted here: https://bitbucket.org/labprentice/splash/src/master/.

## What's new
- Shortwave radiation as input instead of cloudiness.
- Terrain effects on the analytical integrals of the daily energy fluxes.
- Daily infiltration as an analytical integral of the Green-Amp model with corrections for slope.
- Dunne and/or Hortonian runoff generation.
- Analytical solutions for lateral flow and soil water content at any depth (max 2m. at the moment).
- Soil hydrophysical properties estimated by using globally recalibrated pedotransfer functions.
- Maximum water retention in the soil-column computed by equilibrating gravity pushing down and capillarity pulling up.
- Water viscosity effects on the hydraulic conductivity.
- Implementation of empirical formulations (from global studies) to estimate snowfall occurrence and rainfall/snowfall fractions.
- Snowpack balance calculations.

## Installation
To install the development release  of the `rsplash` package please run: 
```r
if(!require(devtools)){install.packages(devtools)}
devtools::install_github( "dsval/rsplash")
```
## Example
This example runs splash with data from one station of the SNOTEL network (https://wcc.sc.egov.usda.gov/nwcc/site?sitenum=361)
*Soil Water content here (SWC (mm)) is defined as the accumulated soil moisture $\theta$ (v/v) from all the measured depth intervals.
```r
library(rsplash)
# load some data
data(Bourne)
# run splash
run1<-splash.point(
	sw_in=Bourne$forcing$sw_in,	# shortwave radiation W/m2
	tc=Bourne$forcing$Ta,		# air temperature C
	pn= Bourne$forcing$P,		# precipitation mm
	lat=Bourne$md$latitude,		# latitude deg
	elev=Bourne$md$elev_m,		# elevation masl
	slop=Bourne$md$slop_250m,	# slope deg
	asp=Bourne$md$asp_250m,		# aspect deg
	soil_data=Bourne$soil, 		# soil data: sand,clay,som in w/w %. Gravel v/v %, bulk density g/cm3, and depth to the bedrock (m)**
	Au=Bourne$md$Aups_250m,		# upslope area m2
	resolution=250.0  			# resolution pixel dem used to get Au
)
#*NOTE: if slop=0.0 (flat surface) the lateral flow is assumed negligible, so: asp,Au and resolution can be ommitted, it won't affect the calculations since all the fluxes are assumed vertical.
#**Soil column thickness

# plot the snow water equivalent
plot(Bourne$forcing$swe,main='SWE (mm)');lines(run1$snow,col=2,lwd=2)
addLegend(legend.loc = "topright", legend.names =c('SWE obs.','SWE sim.'),col=c(1,2),lty=rep(1,2), lwd=rep(2,2))

#Compare the simulations of soil water content (mm) with the measurements taken up to Bourne$max_depth_sm (0.49 m):
# get the simulated water content in the measured region of the profile
swc<-unSWC(soil_data = Bourne$soil, uns_depth = Bourne$max_depth_sm,wn = run1$wn)

# plot the soil water content up to 0.49 m
dev.new()
plot(Bourne$forcing$sm,main=paste('SWC (mm)','up to',Bourne$max_depth_sm,'m'))
lines(swc[[1]],col=4,lwd=2)
addLegend(legend.loc = "topright", legend.names =c('SWC obs.','SWC sim.'),col=c(1,4),lty=rep(1,2), lwd=rep(2,2))

```

## References

Sandoval, D., Prentice, I. C., et. al., Simple process-led algorithms for simulating habitats (SPLASH v.2.0): Robust calculations of water and energy fluxes (in progress)

Davis, T.W., Prentice, I.C., Stocker, B.D., Thomas, R.T., Whitley, R.J., Wang, H., Evans, B.J., Gallego-Sala, A. V., Sykes, M.T., Cramer, W., 2017. Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of radiation, evapotranspiration and plant-available moisture. Geosci. Model Dev. 10, 689–708. doi:10.5194/gmd-10-689-2017
