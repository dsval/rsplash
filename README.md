## Overview

`rsplash` is the R implementation of the 

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
	slop=Bourne$md$slop_250m,	# slope deg
	asp=Bourne$md$asp_250m,		# aspect deg
	soil_data=Bourne$soil, 		# soil data: sand,clay,som,garvel, all in %: bulk density g/cm3
	Au=Bourne$md$Aups_250m,		# upslope area m2
	resolution=250.0  			# resolution pixel dem used to get Au
)
# plot the snow water equivalent
plot(Bourne$forcing[,4],main='SWE (mm)');lines(run1[,5],col=2,lwd=2)

#To compare the simulations of soil water content (mm) with measurements usually taken up to a certain depth:

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
