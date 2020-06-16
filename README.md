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
- Implementation of empirical formulations (from global studies) for snowfall ocurrence and rainfall/snowfall fraction.
- Snowpack balance calculations.

## Installation
To install the development release  of the `rsplash` package please run: 
```r
devtools::install_github( "dsval/rsplash")
```
## Example
This example runs splash with data from one station of the SNOTEL network (https://wcc.sc.egov.usda.gov/nwcc/site?sitenum=361)
*Soil Water content here (SWC) is defined as the accumulated soil moisture (volumetric) from all the depth intervals.
```r
library(rsplash)
# load some data
data(Bourne)
# run splash
run1<-splash.point(
	sw_in=Bourne$forcing[,3],	# shortwave radiation W/m2
	tc=Bourne$forcing[,2],		# air temperature C
	pn= Bourne$forcing[,1],		# precipitation mm
	lat=Bourne$md$latitude,		# latitude deg
	elev=Bourne$md$elev_m,		# elevation masl
	slop=Bourne$md$slop_250m,	# slope* deg 
	asp=Bourne$md$asp_250m,		# aspect deg
	soil_data=Bourne$soil, 		# soil data: sand,clay,som,gravel, all in %; bulk density g/cm3 and depth** (m)
	Au=Bourne$md$Aups_250m,		# upslope area m2
	resolution=250.0  			# resolution of the dem used to get Au
)
#*NOTE: if slop=0.0 (flat surface) the lateral flow is assumed negligible, so: asp,Au and resolution can be ommitted, it won't affect the calculations since all the fluxes are assumed vertical.
#**Soil column thickness

# plot the snow water equivalent
plot(Bourne$forcing[,4],main='SWE (mm)');lines(run1[,5],col=2,lwd=2)

#Compare the simulations of soil water content (mm) with measurements taken up to Bourne$max_depth_sm:

# get the water content in the measured region of the profile
swc<-unSWC(soil_data = Bourne$soil, uns_depth = Bourne$max_depth_sm,wn = run1[,1])
# plot soil water content
dev.new()
plot(Bourne$forcing[,5],main=paste('SWC (mm)','up to',Bourne$max_depth_sm,'m'))
lines(swc[[1]],col=4)

```

## References

Sandoval, D., Prentice, I. C., et. al., Simple process-led algorithms for simulating habitats (SPLASH v.2.0): Robust calculations of water and energy fluxes (in progress)

Davis, T.W., Prentice, I.C., Stocker, B.D., Thomas, R.T., Whitley, R.J., Wang, H., Evans, B.J., Gallego-Sala, A. V., Sykes, M.T., Cramer, W., 2017. Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of radiation, evapotranspiration and plant-available moisture. Geosci. Model Dev. 10, 689–708. doi:10.5194/gmd-10-689-2017
