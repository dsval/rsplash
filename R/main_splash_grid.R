# main splash point
#
# VERSION: 2.0
# LAST UPDATED: 2012-02-19
#
# ~~~~~~~~
# license:
# ~~~~~~~~
# Copyright (C) 2016 Prentice Lab
#
# This file is part of the SPLASH model.
#
# SPLASH is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 2.1 of the License, or
# (at your option) any later version.
#
# SPLASH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SPLASH.  If not, see <http://www.gnu.org/licenses/>.
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led
# algorithms for simulating habitats (SPLASH): Robust indices of radiation
# evapo-transpiration and plant-available moisture, Geoscientific Model
# Development, 2016 (in progress)
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script runs the SPLASH model for one year.
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# - updated monthly and daily results & process [15.01.13]
# - updated plots of results [15.01.16]
# - added example data CSV file [15.01.16]
# - fixed Cramer-Prentice alpha definition [15.01.16]
# - updated monthly results plot [15.01.22]
# - added write out for daily/monthly results [15.01.27]
# - added example of yearly looping [15.01.27]
# - updated year extraction from filename in yearly loop example [15.01.30]
# - fixed plots for Figs. 3 & 4 of manuscript [15.11.23]
# - fixed directory paths [16.02.17]
# - working with any grid datasets [16.02.18]
# - spin-up pixels in parallel [16.02.18]
# - call topmodel R package to calculate topo idx, streams and upslope area  [16.02.18]
# - extract latitude from the data  [16.02.18]
# - calculates slope and aspect from "raster" R package  [16.02.18]
# - calculates streamflow from runoff and baseflow grids

require(Rcpp)
Rcpp::loadModule("splash_module", TRUE)

#### DEFINE FUNCTIONS ########################################################
# Hydrophysics V2
# ************************************************************************
# Name:     soil_hydro
# Input:    - float, sand, (percent)
#           - float, clay, (percent)
#           - float, OM Organic Matter (percent)
#           - float,fgravel, (percent-volumetric)
# Output:   list:
#           - float, FC, (volumetric fraction)
#           - float, WP (volumetric fraction)
#           - float,SAT, (volumetric fraction)
#           - float, AWC (volumetric fraction)
#           - float,Ksat, Saturate hydraulic conductivity/infiltration capacity(mm/hr)
#           - float, A (Coefficient)
#           - float, B (Clapp and Hornberger (1978) pore-size distribution index)
# Features: calculate some soil hydrophysic characteristics
# Ref:      Saxton, K.E., Rawls, W.J., 2006. Soil Water Characteristic Estimates 
#           by Texture and Organic Matter for Hydrologic Solutions. 
#           Soil Sci. Soc. Am. J. 70, 1569. doi:10.2136/sssaj2005.0117
# ************************************************************************

soil_hydro<-function(sand, clay, OM, fgravel=0,bd=NA, ...) {
	
	results<-list()
	# testing
	# sand<-60
	# clay<-30
	# silt<-100-sand-clay
	# OM<-10
	# end test
	# get fractions 
	sand<-sand/100
	clay<-clay/100
	OM<-OM/100
	fgravel<-fgravel/100
	
	depth<-30
	dp<-1/((OM/1.3)+((1-OM)/2.65))
	
	if(is.na(bd)){
		bd<-(1.5 + (dp-1.5-1.10*(1 - clay))*(1-exp(-0.022*depth)))/(1+6.27*OM)
	} 
	
	
	sat<-1-(bd/dp)
	
	# fc<-(sat/bd)*(0.565 + (0.991 - 0.565)*clay^0.5)*exp(-(0.103*sand - 0.785* OM)/(sat/bd))
	fc<-(sat/bd)*(0.3366685 + (1.417544 - 0.3366685)*clay^0.5)*exp(-(0.03320495*sand - 0.2755312* OM)/(sat/bd))
	# fc<-
	fc[fc<0]<-0.1
	fc[fc>1]<-1
	
	wp<- fc*(0.1437904 + (0.8398534 - 0.1437904)*clay^0.5) 
	
	L_10_Ksat<- -2.793574+3.12048*log10(dp-bd)+4.358185*sand
	ksat<-10^L_10_Ksat
	# to mm/h
	ksat<-ksat*10
	
	moist_fvol33init<-0.278*sand+0.034*clay+0.022*OM-0.018*(sand*OM)-0.027*(clay*OM)-0.584*(sand*clay)+0.078
	moist_fvol33<-moist_fvol33init+(0.636*moist_fvol33init-0.107)
	
	# get parameters for BC eqn form SAxton 2006
	coef_B<-(log(1500)-log(33))/(log(fc)-log(wp))
	coef_A<-exp(log(33)+coef_B*log(fc))
	coef_lambda<-1/coef_B
	# Ksat<-1930*(SAT_fvol-FC_fvol)^(3-coef_lambda)
	
	bub_init<--21.6*sand-27.93*clay-81.97*moist_fvol33+71.12*(sand*moist_fvol33)+8.29*(clay*moist_fvol33)+14.05*(sand*clay)+27.16
	bubbling_p<-bub_init+(0.02*bub_init^2-0.113*bub_init-0.7)
	# 101.97162129779 converts from KPa to mmH2O
	bubbling_p<-bubbling_p*-101.97162129779
	# error in empirical fitting, not possible matric potential positive
	# bubbling_p<-ifelse(bubbling_p>0,bubbling_p*-1,bubbling_p)
	bubbling_p[bubbling_p>0]<-bubbling_p[bubbling_p>0]*-1
	# residual water content for BC eqn, Rawls, 1985
	sand<-sand*100
	clay<-clay*100
	silt<-100-sand-clay
	OM<-OM*100
	
	# Ksat<-10*2.54*10^(-0.6+0.012*sand-0.0064*clay)
	RES<--0.018+0.0009*sand+0.005*clay+0.029*sat -0.0002*clay^2-0.001*sand*sat-0.0002*clay^2*sat^2+0.0003*clay^2*sat -0.002*sat^2*clay
	
	# parameters for van Genutchen eqn
	topsoil<-1
	
	alpha<-exp(-14.96 + 0.03135*clay + 0.0351*silt + 0.646*OM +15.29*dp - 0.192*topsoil -4.671*dp^2- 0.000781*clay^2 - 0.00687*OM^2 + 0.0449/OM + 0.0663*log(silt) + 0.1482*log(OM) - 0.04546*dp *silt - 0.4852*dp*OM + 0.00673*topsoil*clay)
	
	n<-1.0+exp(-25.23 - 0.02195*clay + 0.0074*silt - 0.1940*OM + 45.5*dp - 7.24*dp^2 +0.0003658*clay^2 + 0.002885*OM^2 -12.81/dp - 0.1524/silt - 0.01958/OM - 0.2876*log(silt) - 0.0709*log(OM) -44.6*log(dp) - 0.02264*dp*clay + 0.0896*dp*OM +0.00718*topsoil*clay)
	
	m<-1-(1/n)
	
	
	results$SAT<-sat*(1-fgravel)
	results$FC<-fc*(1-fgravel)
	results$WP<-wp*(1-fgravel)
	results$bd<-bd
	results$AWC<-(fc-wp)
	results$Ksat<-ksat
	results$A<-coef_A
	results$B<-coef_B
	results$RES<-RES*(1-fgravel)
	results$bubbling_p<-bubbling_p
	results$VG_alpha<-alpha
	results$VG_n<-n
	
	# results$VG_m<-m
	
	return(results)
}
# ************************************************************************
# Name:     julian_day
# Inputs:   - double, year (y)
#           - double, month (m)
#           - double, day of month (i)
# Returns:  double, Julian day
# Features: This function converts a date in the Gregorian calendar
#           to a Julian day number (i.e., a method of consecutative
#           numbering of days---does not have anything to do with
#           the Julian calendar!)
#           * valid for dates after -4712 January 1 (i.e., jde >= 0)
# Ref:      Eq. 7.1 J. Meeus (1991), Chapter 7 "Julian Day", Astronomical
#             Algorithms
# ************************************************************************
julian_day <- function(y, m, i) {
	if (m <= 2) {
		y <- y - 1
		m <- m + 12
	}
	a <- floor(y/100)
	b <- 2 - a + floor(a/4)
	
	jde <- floor(365.25*(y + 4716)) + floor(30.6001*(m + 1)) + i + b - 1524.5
	return(jde)
}
# ************************************************************************
# Name:     dsin
# Inputs:   double (d), angle in degrees
# Returns:  double, sine of angle
# Features: This function calculates the sine of an angle (d) given
#           in degrees.
# Depends:  pir
# ************************************************************************
dsin <- function(d) {
	pir <- pi/180       # pi in radians
	sin(d*pir)
}
# ************************************************************************
# Name:     cum.interp, avg.interp
# Inputs:   
#               x.month ..... double, variable monthly value
#               y ....... year
# Returns:  double, daily values
# Features: Interpolate monthly values to daily values
# Depends:  julian day
# ************************************************************************
cum.interp<-function(x.months,y){
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	ndaysmonth<-rep(NA,12)
	for(i in 1: 12){ndaysmonth[i]<-julian_day(y,i+1,1)-julian_day(y,i,1)}
	x.days<-x.months/ndaysmonth
	ind.month<-seq(as.Date(paste(y,1,sep="-"),format="%Y-%j"),as.Date(paste(y,ny,sep="-"),format="%Y-%j"), by="month")
	ind.day<-seq(as.Date(paste(y,1,sep="-"),format="%Y-%j"),as.Date(paste(y,ny,sep="-"),format="%Y-%j"), by="day")
	if (sum(!is.na(x.months)) < 2) {return(rep(NA, ny))} 
	else {approx(ind.month, x.days, ind.day, method = "linear", rule=2)$y}
}

avg.interp<-function(x.months,y){
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	ind.month<-seq(as.Date(paste(y,1,sep="-"),format="%Y-%j"),as.Date(paste(y,ny,sep="-"),format="%Y-%j"), by="month")
	ind.day<-seq(as.Date(paste(y,1,sep="-"),format="%Y-%j"),as.Date(paste(y,ny,sep="-"),format="%Y-%j"), by="day")
	if (sum(!is.na(x.months)) < 2) {return(rep(NA, ny))} 
	else {approx(ind.month, x.months, ind.day, method = "linear", rule=2)$y}
}
# ************************************************************************
# Name:     frain_func
# Inputs:   
#               tc ..... double, air daily temperature C
#			 Tt... double, parameter threshold temperature,where 50% of precipitation falls as rain
#			 Tr... double, TR is range of temperatures where both rainfall and snowfall can occur, in °C, typically around 13 °C
#               y ....... year
# Returns:  fraction of the precipitaion falling as rain
# Features: calculates the fraction of the precipitaion falling as rain, accounting for monthly variations of the parameters 
# Depends:  julian day
# Ref:      Kienzle, 2008
# ************************************************************************
frain_func<-function(tc,Tt,Tr,y,ny=NULL){
	
	if(length(tc)>1){
		ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
		ndaysmonth<-rep(NA,12)
		for(i in 1: 12){ndaysmonth[i]<-julian_day(y,i+1,1)-julian_day(y,i,1)}
		m_ind<-rep(1,ndaysmonth[1])
		for(i in 2:12){m_ind<-c(m_ind,rep(i,ndaysmonth[i]))}	
	}else{
		m_ind<-as.Date(paste(y,ny,sep="-"),format="%Y-%j")
		m_ind<-format(m_ind,"%m")
		m_ind<-as.numeric(m_ind)
	}	
	
	Ttm<-Tt+(Tt*dsin((m_ind+2)/1.91))
	Trm<-Tr*(0.55+dsin(m_ind+4))*0.6
	
	
	frain<-ifelse(tc<=Ttm,5*((tc-Ttm)/(1.4*Trm))^3+6.76*((tc-Ttm)/(1.4*Trm))^2+3.19*((tc-Ttm)/(1.4*Trm))+0.5, 5*((tc-Ttm)/(1.4*Trm))^3-6.76*((tc-Ttm)/(1.4*Trm))^2+3.19*((tc-Ttm)/(1.4*Trm))+0.5)
	frain[frain<0]<-0
	frain[frain>1]<-1
	result<-list(frain,Ttm)
	return(result)
}
# ************************************************************************
snowfall_prob<-function(tc,lat,elev){
	
	# ************************************************************************
	# Name:     snowfall_prob
	# Input:    - float, tc, (deg C)
	
	# Output:   
	#           - float, snowfall occurrence probability
	
	# Features: calculate the probability of snowfall occurrence probability, if >0.5, snowfall will occur
	# Ref:      Jennings, K.S., Winchell, T.S., Livneh, B., Molotch, N.P., 2018. Spatial variation of the 
	# 		  rain-snow temperature threshold across the Northern Hemisphere. Nat. Commun. 9, 1–9. doi:10.1038/s41467-018-03629-7
	# ************************************************************************
	p_snow<-1/(1+exp(-0.5827+1.319*tc-elev*4.18E-4-abs(lat)*1.140E-2))
	return(p_snow)
}
# ************************************************************************
# ************************************************************************
# ************************************************************************
# ************************************************************************
# Name:     spin_up
# Inputs:   - list, meteorological data (mdat)
#               $num_lines ..... double, length of meteorol. variable lists
#               $lat_deg ....... double latitude (degrees)
#               $elev_m ......... double, elevation (m)
#               $year .......... double, year
#               $sw_in ............ list, fraction of sunshine hours
#               $tair .......... list, mean daily air temperature (deg. C)
#               $pn ............ list, precipitation (mm/d)
#           - list, daily totals (dtot)
#               $wm ............ list, daily soil moisture (mm)
# Returns:  list, daily totals
# Features: Updates the soil moisture in daily totals until equilibrium
# Depends:  quick_run
# ************************************************************************
rspin_up <-function(lat,elev, sw_in, tc, pn, slop,asp, y,soil_data, Au,resolution) {
	# get number of days in the year y
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	# interpolate monthly to daily
	if(length(sw_in)==12){sw_in<-avg.interp(sw_in,y)}
	if(length(tc)==12){tc<-avg.interp(tc,y)}
	if(length(pn)==12){pn<-cum.interp(pn,y)}
	# get soil hydrophysical characteristics
	soil_info<-soil_hydro(sand=soil_data[1],clay=soil_data[2],OM=soil_data[3],fgravel =soil_data[4] ,bd = soil_data[5])
	depth <- soil_data[6]
	SAT<-soil_info$SAT*(1-soil_data[4]/100)*depth*1000
	WP<-soil_info$WP*(1-soil_data[4]/100)*depth*1000
	FC<-soil_info$FC*(1-soil_data[4]/100)*depth*1000
	RES<-soil_info$RES*(1-soil_data[4]/100)*depth*1000
	lambda<-1/soil_info$B
	bub_press<-soil_info$bubbling_p
	soil_info<-c(SAT,WP,FC,soil_info$Ksat,lambda,depth,bub_press,RES,Au,resolution^2)
	# define snowfall occurrence:
	# 1. get snowfall probability of occurrence
	p_snow<-snowfall_prob(tc,lat,elev)
	# 2. get the treshold for snowfall occurrence
	Tt<-max(tc[p_snow>=0.5])
	# 3. get the fraction of precipitation falling as rain
	f_rain<-ifelse(p_snow>=0.5,frain_func(tc,Tt,13.3,y)[[1]],1)
	# define snowfall and rainfall:
	snowfall<-pn*(1-f_rain)
	pn<-pn*f_rain
	# initialize c++ program
	my_splash = new(SPLASH, lat, elev)
	# run spin up
	result<-my_splash$spin_up(as.integer(ny), as.integer(y), as.numeric(sw_in), as.numeric(tc),as.numeric(pn),slop,asp,as.numeric(snowfall),soil_info)
	
	return(result)
}
# ************************************************************************
# Name:     run_one_year
# Inputs:   - double, latitude, deg (lat)
#           - double, elevation, m (elev)
#           - double, day of year (n)
#           - double, year (y)
#           - double, daily soil moisture content, mm (wn)
#           - double, daily fraction of bright sunshine (sw_in)
#           - double, daily air temperature, deg C (tc)
#           - double, daily precipitation, mm (pn)
# Returns:  list
#             $ho - daily solar irradiation, J/m2
#             $hn - daily net radiation, J/m2
#             $ppfd - daily PPFD, mol/m2
#             $cond - daily condensation water, mm
#             $eet - daily equilibrium ET, mm
#             $pet - daily potential ET, mm
#             $aet - daily actual ET, mm
#             $wn - daily soil moisture, mm
#             $ro - daily runoff, mm
# Features: Runs SPLASH at a single location/pixel for one year.
# Depends:  run_one_day
# ************************************************************************
run_one_year <- function(lat,elev,slop,asp,sw_in, tc, pn, wn, y, snow,soil_data,Au,resolution) {
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	if(length(sw_in)==12){sw_in<-avg.interp(sw_in,y)}
	if(length(tc)==12){tc<-avg.interp(tc,y)}
	if(length(pn)==12){pn<-cum.interp(pn,y)}
	if(length(wn)==12){wn<-avg.interp(wn,y)}
	
	if(length(wn)>length(pn)){
		wn<-wn[1:length(pn)]
		snow<-snow[1:length(pn)]
	}else if(length(wn)<length(pn)){
		wn<-c(wn,wn[length(wn)])
		snow<-c(snow,snow[length(snow)])
	}
	
	
	# get soil hydrophysical characteristics
	soil_info<-soil_hydro(sand=soil_data[1],clay=soil_data[2],OM=soil_data[3],fgravel =soil_data[4] ,bd = soil_data[5])
	depth <- soil_data[6]
	SAT<-soil_info$SAT*(1-soil_data[4]/100)*depth*1000
	WP<-soil_info$WP*(1-soil_data[4]/100)*depth*1000
	FC<-soil_info$FC*(1-soil_data[4]/100)*depth*1000
	RES<-soil_info$RES*(1-soil_data[4]/100)*depth*1000
	lambda<-1/soil_info$B
	bub_press<-soil_info$bubbling_p
	soil_info<-c(SAT,WP,FC,soil_info$Ksat,lambda,depth,bub_press,RES,Au,resolution^2)
	
	# define snowfall occurrence:
	# 1. get snowfall probability of occurrence
	p_snow<-snowfall_prob(tc,lat,elev)
	# 3. get the treshold for snowfall occurrence
	if(length(tc[p_snow>=0.5])>=1){
		Tt<-max(tc[p_snow>=0.5])	
	}else{
		Tt<-0
	}
	
	# 2. get the fraction of precipitation falling as rain
	f_rain<-ifelse(p_snow>=0.5,frain_func(tc,Tt,13.3,y)[[1]],1)
	# define snowfall and rainfall:
	snowfall<-pn*(1-f_rain)
	pn<-pn*f_rain
	# initialize c++ program
	my_splash = new(SPLASH, lat, elev)
	# run spin up
	
	daily_totals<-my_splash$run_one_year(as.integer(ny), as.integer(y),as.numeric(sw_in),as.numeric(tc),as.numeric(pn),as.numeric(wn),slop,asp,as.numeric(snow),as.numeric(snowfall),soil_info)
	return(daily_totals)
}




splash.point<-function(sw_in, tc, pn, lat,elev,slop,asp,soil_data,Au,resolution){
	
	require(xts)
	# Extract time info from data
	# year
	y<-as.numeric(unique(format(time(pn),'%Y')))
	# ndays in the year
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	# time index
	ztime<-time(pn)
	# time frequency
	time.freq<-abs(as.numeric(ztime[1]-ztime[2], units = "days"))
	
	
	
	if (time.freq<2){		
		
		if (length(y)==1){
			initial<-rspin_up(lat,elev, sw_in, tc, pn, slop,asp, y[1],soil_data,Au,resolution)
			
			result<-run_one_year(lat,elev,slop,asp,sw_in, tc, pn,initial$sm, y[1], initial$snow,soil_data,Au,resolution)
			# result<-xts(result,ztime)
			result<-do.call(cbind,result)
		}
		else if(length(y)>1){
			
			end<-cumsum(ny)
			start<-end+1
			result<-list()
			sw_av<-tapply(sw_in,format(time(sw_in),"%j"),mean)
			tc_av<-tapply(tc,format(time(sw_in),"%j"),mean)
			pn_av<-tapply(pn,format(time(sw_in),"%j"),mean)
			# initial<-rspin_up(lat,elev, sw_in[1:ny[1]], tc[1:ny[1]], pn[1:ny[1]], slop,asp, y[1],soil_data,Au,resolution)
			initial<-rspin_up(lat,elev, sw_av, tc_av, pn_av, slop,asp, y[1],soil_data,Au,resolution)
			result[[1]]<-run_one_year(lat,elev,slop,asp,sw_in[1:ny[1]], tc[1:ny[1]],  pn[1:ny[1]],initial$sm, y[1], initial$snow,
				soil_data,Au,resolution)
			
			for (i in 2:length(y)){
				
				stidx<-i-1
				# correct for leap years	
				result[[i]]<-run_one_year(lat,elev,slop,asp, sw_in[start[stidx]:end[i]], tc[start[stidx]:end[i]], pn[start[stidx]:end[i]],
					result[[stidx]]$wn,y[i],result[[stidx]]$snow,soil_data,Au,resolution)
			}
			result<-lapply(result,FUN=as.data.frame)
			result<-do.call(rbind,result)
						
		}
		
	}
	# order results as time series
	
	
	result<-xts(result,ztime)		
	return(result)
}


# main splash grid
spinup.grid<-function(sw_in, tc, pn, elev,lat, terraines,soil, y, resolution,  Au ,inmem=FALSE,outdir=getwd()){
	###############################################################################################
	# 01. create array for equilibrium soil moisture wneq and snow
	###############################################################################################
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	wneq<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(wneq)<-extent(elev)
	snoweq<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(snoweq)<-extent(elev)	
	###############################################################################################
	# 00. ifferent strategies for large rasters
	###############################################################################################
	if (inmem==TRUE){
		###############################################################################################
		# 00. Preprocess send everything to the Ram
		###############################################################################################
		sw_in<-getValues(sw_in)
		pn<-getValues(pn)
		elev<-getValues(elev)
		tc<-getValues(tc)
		lat<-getValues(lat)
		terraines<-getValues(terraines)
		soil<-getValues(soil)
		Au<-getValues(Au)
		###############################################################################################
		# 02. Loop spin up through pixels in parallel
		###############################################################################################
		niter<-ncell(elev)
		# library(doSNOW)
		# cl <- makeCluster(7)
		# registerDoSNOW(cl)
		pb <- txtProgressBar(max = niter, style = 3)
		progress <- function(n) setTxtProgressBar(pb, n)
		opts <- list(progress = progress)
		
		wneqmat<-foreach (i = icount(niter),.packages = c("raster","Rcpp","rsplash"),.combine=rbind,.multicombine=TRUE,
						.maxcombine=niter,.inorder=TRUE,.options.snow = opts) %dopar% {
			
			rspin_up(lat[i],elev[i],as.numeric(sw_in[i,]),as.numeric(tc[i,]),as.numeric(pn[i,]),
				as.numeric(terraines[i,1]),as.numeric(terraines[i,2]),y,as.numeric(soil[i,]), Au[i],resolution)
			
		}
		# stopCluster(cl)
		
		###############################################################################################
		# 03. Allocate values to raster 
		###############################################################################################
		
		wneq<-setValues(wneq,do.call(rbind,wneqmat[,1]))
		snoweq<-setValues(snoweq,do.call(rbind,wneqmat[,2]))
		rm(wneqmat)
		gc()	
		
	}else{
		###############################################################################################
		# 02. Loop spin up through pixels in parallel
		###############################################################################################
		# library(doSNOW)
		# cl <- makeCluster(6)
		# registerDoSNOW(cl)
		setwd(dirname(rasterTmpFile()))
		wneq<-writeStart(wneq,filename="wneq.grd",overwrite=TRUE)
		snoweq<-writeStart(snoweq,filename="snoweq.grd",overwrite=TRUE)
		pb <- txtProgressBar(min=1,max = nrow(elev), style = 3)
		for(nr in 1:nrow(elev)){
			setTxtProgressBar(pb,nr)
			swrow<-getValues(sw_in,nr)
			tcrow<-getValues(tc,nr)
			pnrow<-getValues(pn,nr)
			latrow<-getValues(lat,nr)
			elevrow<-getValues(elev,nr)
			terrainesrow<-getValues(terraines,nr)
			soilrow<-getValues(soil,nr)
			Aurow<-getValues(Au,nr)
			###############################################################################################
			# 00. Loop spin up through rows in parallel and write values on the go
			###############################################################################################	
			niter<-length(elevrow)
						
			wneqmat<-foreach (i = icount(niter),.packages = c("raster","Rcpp","rsplash"),.combine=rbind,.multicombine=TRUE,.maxcombine=niter,.inorder=TRUE) %dopar% {
				
				rspin_up(latrow[i],elevrow[i],as.numeric(swrow[i,]),as.numeric(tcrow[i,]),as.numeric(pnrow[i,]),
					as.numeric(terrainesrow[i,1]),as.numeric(terrainesrow[i,2]),y,as.numeric(soilrow[i,]), Aurow[i],resolution)
				
			}
			wneq <- writeValues(wneq,do.call(rbind,wneqmat[,1]), nr)
			snoweq <- writeValues(snoweq, do.call(rbind,wneqmat[,2]), nr)
			
		}
		close(pb)
		wneq <- writeStop(wneq)
		snoweq <- writeStop(snoweq)	
		# stopCluster(cl)
		rm(wneqmat)
		gc()		
	}

	return(list(wneq=wneq,snoweq=snoweq))
		
}


run_one_year.grid<-function(sw_in, tc, pn,wn,snow ,elev,lat, terraines,soil, y, resolution,  Au,inmem=FALSE){
	###############################################################################################
	# 00. create array for results, fluxes: mm/day, storages (wn, snow): mm
	###############################################################################################
	
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	# actual soil moisture
	sm<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(sm)<-extent(elev)
	# runoff
	ro<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(ro)<-extent(elev)
	# Potential evapotranspiration
	pet<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(pet)<-extent(elev)
	# Actual evapotranspiration
	aet<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(aet)<-extent(elev)
	# Snow water equivalent
	swe<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(swe)<-extent(elev)
	# Condensation
	cond<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(cond)<-extent(elev)
	# baseflow
	bflow<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(bflow)<-extent(elev)
	if (inmem==TRUE){
		###############################################################################################
		# 01. Preprocess: get vectors and matrices, send everything to the RAM
		###############################################################################################
		sw_in<-getValues(sw_in)
		pn<-getValues(pn)
		wn<-getValues(wn)
		snow<-getValues(snow)
		elev<-getValues(elev)
		tc<-getValues(tc)
		lat<-getValues(lat)
		terraines<-getValues(terraines)
		soil<-getValues(soil)
		Au<-getValues(Au)
		
		###############################################################################################
		# 02. Loop run_one_year through pixels in parallel
		###############################################################################################
		niter<-ncell(elev)
		# library(doSNOW)
		# cl <- makeCluster(7)
		# registerDoSNOW(cl)
		# create function that will separate out foreach output into list of three lists
		
		yearlist<-foreach (i = icount(niter),.packages = c("raster","Rcpp","rsplash"),.combine=rbind,.multicombine=TRUE,.maxcombine=niter,.inorder=TRUE) %dopar% {
			
			run_one_year(lat[i],elev[i],as.numeric(terraines[i,1]),as.numeric(terraines[i,2]),as.numeric(sw_in[i,]),
						as.numeric(tc[i,]),as.numeric(pn[i,]),as.numeric(wn[i,]),y, as.numeric(snow[i,]), as.numeric(soil[i,]), Au[i],resolution)
		}
		# stopCluster(cl)
		
		###############################################################################################
		# 03. Allocate values to raster
		###############################################################################################
		sm<-setValues(sm,do.call(rbind,yearlist[,1]))
		# runoff
		ro<-setValues(ro,do.call(rbind,yearlist[,2]))
		# Potential evapotranspiration
		pet<-setValues(pet,do.call(rbind,yearlist[,3]))
		# Actual evapotranspiration
		aet<-setValues(aet,do.call(rbind,yearlist[,4]))
		# Snow water equivalent
		swe<-setValues(swe,do.call(rbind,yearlist[,5]))
		# Condensation
		cond<-setValues(cond,do.call(rbind,yearlist[,6]))
		# baseflow
		bflow<-setValues(bflow,do.call(rbind,yearlist[,7]))
		rm(yearlist)
		gc()
	}else{
		###############################################################################################
		# 02. Loop spin up through pixels in parallel
		###############################################################################################
		# library(doSNOW)
		# cl <- makeCluster(6)
		# registerDoSNOW(cl)
		setwd(dirname(rasterTmpFile()))
		sm<-writeStart(sm,filename="sm.grd",overwrite=TRUE)
		ro<-writeStart(ro,filename="ro.grd",overwrite=TRUE)
		pet<-writeStart(pet,filename="pet.grd",overwrite=TRUE)
		aet<-writeStart(aet,filename="aet.grd",overwrite=TRUE)
		swe<-writeStart(swe,filename="swe.grd",overwrite=TRUE)
		cond<-writeStart(cond,filename="cond.grd",overwrite=TRUE)
		bflow<-writeStart(bflow,filename="bflow.grd",overwrite=TRUE)
		
		for(nr in 1:nrow(elev)){
			
			swrow<-getValues(sw_in,nr)
			tcrow<-getValues(tc,nr)
			pnrow<-getValues(pn,nr)
			latrow<-getValues(lat,nr)
			elevrow<-getValues(elev,nr)
			wnrow<-getValues(wn,nr)
			snowrow<-getValues(snow,nr)
			terrainesrow<-getValues(terraines,nr)
			soilrow<-getValues(soil,nr)
			Aurow<-getValues(Au,nr)
			###############################################################################################
			# 00. Loop spin up through rows in parallel and write values on the go
			###############################################################################################	
			niter<-length(elevrow)
			
			yearlist<-foreach (i = icount(niter),.packages = c("raster","Rcpp","rsplash"),.combine=rbind,.multicombine=TRUE,.maxcombine=niter,.inorder=TRUE) %dopar% {
				
				run_one_year(latrow[i],elevrow[i],as.numeric(terrainesrow[i,1]),as.numeric(terrainesrow[i,2]),as.numeric(swrow[i,]),
					as.numeric(tcrow[i,]),as.numeric(pnrow[i,]),as.numeric(wnrow[i,]),y, as.numeric(snowrow[i,]), as.numeric(soilrow[i,]), Aurow[i],resolution)
				
			}
			sm <- writeValues(sm,do.call(rbind,yearlist[,1]), nr)
			ro <- writeValues(ro, do.call(rbind,yearlist[,2]), nr)
			pet <- writeValues(pet, do.call(rbind,yearlist[,3]), nr)
			aet <- writeValues(aet, do.call(rbind,yearlist[,4]), nr)
			swe <- writeValues(swe, do.call(rbind,yearlist[,5]), nr)
			cond <- writeValues(cond, do.call(rbind,yearlist[,6]), nr)
			bflow <- writeValues(bflow, do.call(rbind,yearlist[,7]), nr)
			
		}
		
		sm <- writeStop(sm)
		ro <- writeStop(ro)
		pet <- writeStop(pet)
		aet <- writeStop(aet)
		swe <- writeStop(swe)
		cond <- writeStop(cond)
		bflow <- writeStop(bflow)	
		# stopCluster(cl)
		
		rm(yearlist)
		gc()
		
	}
	
	
	return(list(wn=stack(sm),ro=stack(ro),pet=stack(pet),aet=stack(aet),snow=stack(swe),cond=stack(cond),bflow=stack(bflow)))
	
	
}

splash.grid<-function(sw_in, tc, pn, elev, soil, outdir=getwd(),sim.control=list(par=TRUE, ncores=7,output.mode="monthly",inmem=FALSE), ...){
	#### IMPORT SOURCES ##########################################################
	require(raster)
	require(xts)
	require(doSNOW)
	require(zoo)
	rasterOptions(maxmemory=3e7, tmptime = 24, chunksize = 3e7,todisk = FALSE, overwrite=TRUE, tolerance = 0.5)
	
	###########################################################################
	# 01. Calculate spatial distributed variables
	###########################################################################
	setwd(dirname(rasterTmpFile()))
	# get resolution in m2
	resolution<-sqrt(cellStats(area(elev), stat='mean', na.rm=TRUE))*1000
	# 1.1 calculate upslope area in m2, *all the raster will be saved to the disk by default
	if (ncell(elev)>1.5e7|resolution>=25000){
		Au<-upslope_areav2(elev)
	}else{
		Au<-upslope_area(elev)
	}
	
	# 1.2 get latitudes
	lat<-elev*0
	lat.data<-rasterToPoints(elev)
	lat[!is.na(lat)]<-lat.data[,2]
	rm(lat.data)
	#lat<-writeRaster(lat,filename="lat.grd", overwrite=TRUE)
	
	# 1.3  calculate slope and aspect
	terraines<-terrain(elev, opt=c('slope', 'aspect'), unit='degrees')
	
	###########################################################################
	# 02. Get time info from data
	###########################################################################
	
	y<-as.numeric(unique(format(getZ(pn),'%Y')))
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	ztime<-getZ(pn)
	ztime.months<-seq(as.Date(paste(y[1],1,sep="-"),format="%Y-%j"),as.Date(paste(y[length(y)],ny[length(y)],sep="-"),format="%Y-%j"), by="month")
	ztime.days<-seq(as.Date(paste(y[1],1,sep="-"),format="%Y-%j"),as.Date(paste(y[length(y)],ny[length(y)],sep="-"),format="%Y-%j"), by="day")
	zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31))
	time.freq<-ztime[2]-ztime[1]
	###########################################################################
	# 03. Start the calculations for daily inputs
	###########################################################################
	if (abs(as.numeric(time.freq, units = "days"))<2){		
		
		if (length(y)==1){
			###########################################################################
			# 3.1. Equilibrate
			###########################################################################
			cat("reaching steady state")
			cl <- parallel::makeCluster(sim.control$ncores, ...)
			doSNOW::registerDoSNOW(cl)
			eq<-spinup.grid(sw_in,tc,pn,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir)
			cat(paste("solving","year",y[1]))
			result.all<-run_one_year.grid(sw_in,tc,pn,eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem)
			stopCluster(cl)
		}
		else if(length(y)>1){
			
			end<-cumsum(ny)
			start<-end+1
			result<-list()
			cl <- parallel::makeCluster(sim.control$ncores, ...)
			doSNOW::registerDoSNOW(cl)
			cat("reaching steady state")
			eq<-spinup.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]],elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir)
			cat(paste("solving","year",y[1]))
			result[[1]]<-run_one_year.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]], eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],
										resolution,Au,sim.control$inmem)
			
			###########################################################################
			# 3.2. loop through years
			###########################################################################
			# correct for leap years	inside c++ functions
			pb <- txtProgressBar(min=1,max = length(y), style = 3)
			for (i in 2:length(y)){
				cat(paste("solving","year",y[i]))
				setTxtProgressBar(pb,i)
				stidx<-i-1
				result[[i]]<-run_one_year.grid(sw_in[[start[stidx]:end[i]]], tc[[start[stidx]:end[i]]], pn[[start[stidx]:end[i]]],result[[stidx]]$wn,
					result[[stidx]]$snow,elev,lat,terraines,soil,y[i],resolution,Au,sim.control$inmem)				
			}
			close(pb)
			stopCluster(cl)
			
			
			gc()
			
		}
		
		
	}
	
	###########################################################################
	# 03. Start the calculations for Monthly inputs
	###########################################################################
	
	else if (abs(as.numeric(time.freq, units = "days"))>20){		
		
		if (length(y)==1){
			cat("reaching steady state")
			cl <- parallel::makeCluster(sim.control$ncores, ...)
			doSNOW::registerDoSNOW(cl)
			eq<-spinup.grid(sw_in,tc,pn,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir)
			cat(paste("solving","year",y[1]))
			result.all<-run_one_year.grid(sw_in,tc,pn,eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem)
			stopCluster(cl)
		}
		else if(length(y)>1){
			nm <- rep(12,length(y))
			end<-cumsum(nm)
			start<-end-11
			result<-list()
			cl <- parallel::makeCluster(sim.control$ncores, ...)
			doSNOW::registerDoSNOW(cl)
			cat("reaching steady state")
			eq<-spinup.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]],elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir)
			cat(paste("solving","year",y[1]))
			result[[1]]<-run_one_year.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]], eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],
										resolution,Au,sim.control$inmem)
			
			###########################################################################
			# 3.2. loop through years
			###########################################################################
			# correct for leap years	inside c++ functions
			pb <- txtProgressBar(min=1,max = length(y), style = 3)
			
			for (i in 2:length(y)){
				cat(paste("solving","year",y[i]))
				setTxtProgressBar(pb,i)
				stidx<-i-1
				
				result[[i]]<-run_one_year.grid(sw_in[[start[i]:end[i]]], tc[[start[i]:end[i]]], pn[[start[i]:end[i]]],result[[stidx]]$wn,
					result[[stidx]]$snow,elev,lat,terraines,soil,y[i],resolution,Au,sim.control$inmem)
				
				
			}
			close(pb)
			stopCluster(cl)
			
			gc()
		}
		
	}		
	# 
	###########################################################################
	# 4. Building the raster stacks
	###########################################################################
	if(length(y)>1){
		
		cat("building the grids")
		gc()
		result.all<-list()
		# memory leak???	
		result.all$wn<-result[[1]]$wn
		result.all$wn@layers<-purrr::flatten(lapply(result, function(x) {as.list(x$wn)}))
		result.all$ro<-result[[1]]$ro
		result.all$ro@layers<-purrr::flatten(lapply(result, function(x) {as.list(x$ro)}))
		result.all$pet<-result[[1]]$pet
		result.all$pet@layers<-purrr::flatten(lapply(result, function(x) {as.list(x$pet)}))
		result.all$aet<-result[[1]]$aet
		result.all$aet@layers<-purrr::flatten(lapply(result, function(x) {as.list(x$aet)}))
		result.all$snow<-result[[1]]$snow
		result.all$snow@layers<-purrr::flatten(lapply(result, function(x) {as.list(x$snow)}))
		result.all$cond<-result[[1]]$cond
		result.all$cond@layers<-purrr::flatten(lapply(result, function(x) {as.list(x$cond)}))
		result.all$bflow<-result[[1]]$bflow
		result.all$bflow@layers<-purrr::flatten(lapply(result, function(x) {as.list(x$bflow)}))
		
		rm(result)
	}
	gc()
	
	
	if (sim.control$output.mode=="monthly"){
		cat("aggregating to monthly")
		###########################################################################
		# 5. Aggregate monthly
		###########################################################################
		# define functions to aggregate
		indmonth<-format(ztime.days,'%Y-%m')
		
		month_mean<-function(x){
			
			x<-calc(x=x, fun=function(x)tapply(x,indmonth,FUN= mean))
		}
		month_sum<-function(x){
			x<-calc(x=x, fun=function(x)tapply(x,indmonth,FUN= sum))
		}
		
		result.all[c(1,5)]<-lapply(result.all[c(1,5)], month_mean)
		
		result.all[c(2,3,4,6,7)]<-lapply(result.all[c(2,3,4,6,7)], month_sum)
		# setting ztime
		settime<-function(x){
			x<-raster::setZ(x,ztime.months)
		}
		
		result.all<-lapply(result.all,settime)
		
		gc()
		cat("writing to disk")
		result.all$wn<-writeRaster(result.all$wn,paste0(outdir,"/",y[1],"_",y[length(y)],".","wn",".","nc"),format="CDF",overwrite=TRUE,varname="wn", varunit="mm",
			longname="monthly soil moisture", xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
		result.all$ro<-writeRaster(result.all$ro,paste0(outdir,"/",y[1],"_",y[length(y)],".","ro",".","nc"),format="CDF",overwrite=TRUE,
			varname="ro", varunit="mm", longname="monthly runoff", xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
		result.all$pet<-writeRaster(result.all$pet,paste0(outdir,"/",y[1],"_",y[length(y)],".","pet",".","nc"),format="CDF",overwrite=TRUE,varname="pet", varunit="mm",
			longname="monthly potential ET", xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
		result.all$aet<-writeRaster(result.all$aet,paste0(outdir,"/",y[1],"_",y[length(y)],".","aet",".","nc"),format="CDF",overwrite=TRUE,varname="aet", varunit="mm",								longname="monthly actual ET", xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
		result.all$snow<-writeRaster(result.all$snow,paste0(outdir,"/",y[1],"_",y[length(y)],".","swe",".","nc"),format="CDF",overwrite=TRUE,
			varname="swe", varunit="mm",longname="monthly snow water equivalent", xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
		result.all$cond<-writeRaster(result.all$cond,paste0(outdir,"/",y[1],"_",y[length(y)],".","cond",".","nc"),format="CDF",overwrite=TRUE,varname="cond", 
			varunit="mm", longname="monthly condensation water", xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
		result.all$bflow<-writeRaster(result.all$bflow,paste0(outdir,"/",y[1],"_",y[length(y)],".","bflow",".","nc"),format="CDF",overwrite=TRUE,varname="bflow", 		varunit="mm", longname="monthly baseflow", xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
		
	}else{
		cat("writing to disk")
		result.all$wn<-writeRaster(result.all$wn,paste0(outdir,"/",y[1],"_",y[length(y)],".","wn",".","nc"),format="CDF",overwrite=TRUE,varname="wn", varunit="mm",longname="daily soil moisture", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31)))
		result.all$ro<-writeRaster(result.all$ro,paste0(outdir,"/",y[1],"_",y[length(y)],".","ro",".","nc"),format="CDF",overwrite=TRUE,
			varname="ro", varunit="mm/day", longname="daily runoff", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31)))
		result.all$pet<-writeRaster(result.all$pet,paste0(outdir,"/",y[1],"_",y[length(y)],".","pet",".","nc"),format="CDF",overwrite=TRUE,varname="pet", varunit="mm/day", longname="daily potential ET", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31)))
		result.all$aet<-writeRaster(result.all$aet,paste0(outdir,"/",y[1],"_",y[length(y)],".","aet",".","nc"),format="CDF",overwrite=TRUE,varname="aet", varunit="mm/day",longname="daily actual ET", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31)))
		result.all$snow<-writeRaster(result.all$snow,paste0(outdir,"/",y[1],"_",y[length(y)],".","swe",".","nc"),format="CDF",overwrite=TRUE,
			varname="swe", varunit="mm",longname="daily snow water equivalent", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31)))
		result.all$cond<-writeRaster(result.all$cond,paste0(outdir,"/",y[1],"_",y[length(y)],".","cond",".","nc"),format="CDF",overwrite=TRUE,varname="cond", 
			varunit="mm/day", longname="daily condensation water", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31)))
		result.all$bflow<-writeRaster(result.all$bflow,paste0(outdir,"/",y[1],"_",y[length(y)],".","bflow",".","nc"),format="CDF",overwrite=TRUE,varname="bflow", 
			varunit="mm/day", longname="daily baseflow", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31)))
		
	}
	
	gc()		
	return(result.all)
}
