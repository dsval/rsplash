#' Simple process-led algorithms for simulating habitats (SPLASH v.2.0)
#' 
#' R/C++ implementation of the SPLASH v.2.0 algorithm (Davis et al., 2017; Sandoval et al., in prep.).
#' 
#' @param   sw_in Incoming shortwave solar radiation (W m-2), timeseries object of monthly or daily averages.
#' @param   tc Air temperature (°C), same timestep as sw_in
#' @param   pn Precipitation (mm), same timestep as sw_in
#' @param   elev Elevation (m.a.s.l)
#' @param   slop Terrain feature: slope inclination (°)
#' @param   asp Terrain feature: slope orientation (°), standard clockwise from 0.0° North
#' @param   soil_data Soil data organized as a vector in the way: c(sand(perc),clay(perc),organic matter(perc),coarse-fragments-fraction(perc), bulk density(g cm-3))
#' @return a time series matrix including:
#' \itemize{
#'         \item \eqn{wn}: Soil water content within the first 2m of depth (mm).
#'         \item \eqn{ro}: Runoff (mm d-1).
#'         \item \eqn{pet}: Potential evapotranspiration (mm d-1).
#'         \item \eqn{aet}: Actual evapotranspiration (mm d-1).
#'         \item \eqn{snow}: Snow water equivalent (mm).
#'         \item \eqn{cond}: Condensation (mm d-1).
#'         \item \eqn{bflow}: Lateral flow (mm d-1).
#'         \item \eqn{netr}: Daytime net radiation (MJ d-1).
#' }
#' @import Rcpp 
#' @import xts
#' @keywords splash, evapotranspiration, soil moisture
#' @export
#' @examples
#' splash.point(sw_in=200, tc=15, pn=10, lat=44,elev=1800,slop=10,asp=270,soil_data=c(sand=44,clay=2,OM=6,fgravel=12))
splash.point<-function(sw_in, tc, pn, lat,elev,slop=0,asp=0,soil_data,Au=0,resolution=250){
	###########################################################################
	# 01.Extract time info from the data
	###########################################################################	
	# year
	y<-as.numeric(unique(format(time(pn),'%Y')))
	# ndays in the year
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	# time index
	ztime<-time(pn)
	# time frequency
	time.freq<-abs(as.numeric(ztime[1]-ztime[2], units = "days"))
	
	###########################################################################
	# 02. Start the calculations for daily inputs
	###########################################################################
	
	if (time.freq<2){		
		
		if (length(y)==1){
			initial<-rspin_up(lat,elev, sw_in, tc, pn, slop,asp, y[1],soil_data,Au,resolution)
			result<-run_one_year(lat,elev,slop,asp,sw_in, tc, pn,initial$sm, y[1], initial$snow,soil_data,Au,resolution,initial$qin,initial$tdrain, initial$snwage)
			result<-do.call(cbind,result)
		}
		else if(length(y)>1){
			
			end<-cumsum(ny)
			start<-end+1
			result<-list()
			# average daily for spin-up
			sw_av<-tapply(sw_in,format(time(sw_in),"%j"),mean, na.rm=TRUE)
			tc_av<-tapply(tc,format(time(sw_in),"%j"),mean, na.rm=TRUE)
			pn_av<-tapply(pn,format(time(sw_in),"%j"),mean, na.rm=TRUE)
			
			initial<-rspin_up(lat,elev, sw_av, tc_av, pn_av, slop,asp, y[1],soil_data,Au,resolution)
			result[[1]]<-run_one_year(lat,elev,slop,asp,sw_in[1:ny[1]], tc[1:ny[1]],  pn[1:ny[1]],initial$sm, y[1], initial$snow,
				soil_data,Au,resolution,initial$qin,initial$tdrain, initial$snwage)
			
			for (i in 2:length(y)){
				
				stidx<-i-1
				# correct for leap years inside the c++ code	
				result[[i]]<-run_one_year(lat,elev,slop,asp, sw_in[start[stidx]:end[i]], tc[start[stidx]:end[i]], pn[start[stidx]:end[i]],
					result[[stidx]]$wn,y[i],result[[stidx]]$snow,soil_data,Au,resolution,result[[stidx]]$qin_prev,result[[stidx]]$tdrain,result[[stidx]]$snwage)
			}
			result<-lapply(result,FUN=base::as.data.frame)
			result<-do.call(rbind,result)
			
		}
		result<-xts(result[,1:9],ztime)
	}
	###########################################################################
	# 03. Start calculations if the  inputs are monthly 
	###########################################################################
	if (time.freq>20){		
		ztime.days<-seq(as.Date(paste(y[1],1,sep="-"),format="%Y-%j"),as.Date(paste(y[length(y)],ny[length(y)],sep="-"),format="%Y-%j"), by="day")
		if (length(y)==1){
			initial<-rspin_up(lat,elev, sw_in, tc, pn, slop,asp, y[1],soil_data,Au,resolution)
			
			result<-run_one_year(lat,elev,slop,asp,sw_in, tc, pn,initial$sm, y[1], initial$snow,soil_data,Au,resolution,initial$qin,initial$tdrain, initial$snwage)
			
			result<-do.call(cbind,result)
		}
		else if(length(y)>1){
			nm <- rep(12,length(y))
			end<-cumsum(nm)
			start<-end-11
			result<-list()
			sw_av<-tapply(sw_in,format(time(sw_in),"%m"),mean, na.rm=TRUE)
			tc_av<-tapply(tc,format(time(sw_in),"%m"),mean, na.rm=TRUE)
			pn_av<-tapply(pn,format(time(sw_in),"%m"),mean, na.rm=TRUE)
			initial<-rspin_up(lat,elev, sw_av, tc_av, pn_av, slop,asp, y[1],soil_data,Au,resolution)
			result[[1]]<-run_one_year(lat,elev,slop,asp,sw_in[1:end[1]], tc[1:end[1]],  pn[1:end[1]],initial$sm, y[1], initial$snow,
				soil_data,Au,resolution,initial$qin,initial$tdrain, initial$snwage)
			
			for (i in 2:length(y)){
				
				stidx<-i-1
				# correct for leap year sinside the c++ code		
				result[[i]]<-run_one_year(lat,elev,slop,asp, sw_in[start[i]:end[i]], tc[start[i]:end[i]], pn[start[i]:end[i]],
					result[[stidx]]$wn,y[i],result[[stidx]]$snow,soil_data,Au,resolution,result[[stidx]]$qin_prev,result[[stidx]]$tdrain,result[[stidx]]$snwage)
			}
			result<-lapply(result,FUN=base::as.data.frame)
			result<-do.call(rbind,result)
			
		}
		result<-xts(result[,1:9],ztime.days)
	}
	####################################################################################################
	# 4. compute soil moisture limitations
	####################################################################################################
	#get the wilting point and bucket size in mm
	soil_water<-soil_hydro(sand=soil_data[1],clay=soil_data[2],OM=soil_data[3],fgravel =0.0 ,bd = soil_data[5])
	wp<-(soil_water$WP/2)*soil_data[6]*1000
	KWm<-(soil_water$FC-(soil_water$WP/2))*(soil_data[6]*1000)
	#get relative soil moisture limitation from 0.0 (at WP) to 1.0 (at FC)
	soil_lim<-(result$wn-wp)/KWm
	#adjust the boundaries, wn goes from ~WP to SAT
	soil_lim[soil_lim<0]<-0
	soil_lim[soil_lim>1]<-1
	result$sm_lim<-soil_lim		
	return(result)
}


# require(Rcpp)
Rcpp::loadModule("splash_module", TRUE)


soil_hydro<-function(sand, clay, OM, fgravel=0,bd=NA, ...) {
	# Hydrophysics V2
	# ************************************************************************
	# Name:     soil_hydro
	# Input:    - float, sand, (percent)
	#           - float, clay, (percent)
	#           - float, OM Organic Matter (percent)
	#           - float, fgravel, (percent-volumetric)
	#           - float, bd, bulk density (g/cm3)
	# Output:   list:
	#           - float, FC, (volumetric fraction)
	#           - float, WP (volumetric fraction)
	#           - float,SAT, (volumetric fraction)
	#           - float, AWC (volumetric fraction)
	#           - float,Ksat, Saturate hydraulic conductivity/ min infiltration capacity(mm/hr)
	#           - float, A, B, Coefficients to fit the 33-1500kPa section of the retention curve from Saxton and Rawls, 2006
	#           - float, bubbling_p, air entry pressure
	#           - float, alpha, n, Shape parameters for the van Genuchten equation for water retentin
	# Features: calculate some soil hydrophysic characteristics
	# Ref:      Saxton, K.E., Rawls, W.J., 2006. Soil Water Characteristic Estimates 
	#           by Texture and Organic Matter for Hydrologic Solutions. 
	#           Soil Sci. Soc. Am. J. 70, 1569. doi:10.2136/sssaj2005.0117
	#		  Balland, V., Pollacco, J.A.P., Arp, P.A., 2008. Modeling soil hydraulic properties for 
	#		  a wide range of soil conditions. Ecol. Modell. 219, 300–316. doi:10.1016/j.ecolmodel.2008.07.009
	# ************************************************************************
	results<-list()
	########################################################################################
	# 01. get fractions
	######################################################################################## 
	sand<-sand/100
	clay<-clay/100
	OM<-OM/100
	fgravel<-fgravel/100
	########################################################################################
	# 02. calc bulk density [g/cm3] in case it is not provided, assumming 30 cm depth Balland et al. (2008)
	########################################################################################
	depth<-30
	dp<-1/((OM/1.3)+((1-OM)/2.65))
	
	if(!is.numeric(sand)){
		if(is.null(bd)){
			bd<-(1.5 + (dp-1.5-1.10*(1 - clay))*(1-exp(-0.022*depth)))/(1+6.27*OM)
		} 
	}else{
		bd<-ifelse(is.na(bd),(1.5 + (dp-1.5-1.10*(1 - clay))*(1-exp(-0.022*depth)))/(1+6.27*OM),bd)
		
		}
	########################################################################################
	# 03. calc volumetric water contents at saturation, 33kPa (fc) and 1500kPa (wp) Balland et al. (2008)
	######################################################################################## 
	# volumetric water content at saturation [m3/m3]
	sat<-1-(bd/dp)
	# volumetric water content at 33kPa [m3/m3]
	fc<-(sat/bd)*(0.4760944 + (0.9402962 - 0.4760944)*clay^0.5)*exp(-1*(0.05472678*sand - 0.01* OM)/(sat/bd))
	##errors in the empirical fitting, switch to Saxton and Rawls (2006)
	
	FCinit<--0.251*sand+0.195*clay+0.011*OM+0.006*(sand*OM)-0.027*(clay*OM)+0.452*(sand*clay)+0.299
	FC_fvol<-FCinit+(1.283*FCinit^2-0.374*FCinit-0.015)
	fc[fc<0 | fc>sat]<-FC_fvol[fc<0 | fc>=sat]
	##persistent errors, few pixels in siberia
	#fc[fc<0]<-0.1
	#fc[fc>sat]<-0.9*sat
	# if(!is.numeric(sand)){
	# 	fc[fc>sat]<-0.9*sat[fc>sat]
	# }else{
	# 	fc[fc>sat]<-0.9*sat
	# }
	#fc[fc>sat]<-0.9*sat[fc>sat]
	# volumetric water content at 1500kPa [m3/m3]
	wp<- fc*(0.2018522 + (0.7809203 - 0.2018522)*clay^0.5) 
	
	
	########################################################################################
	# 05. calc shape parameters water retention Brooks and Corey curve Saxton and Rawls (2006)
	######################################################################################## 
	# get parameters for BC eqn form Saxton 2006
	coef_B<-(log(1500)-log(33))/(log(fc)-log(wp))
	coef_A<-exp(log(33)+coef_B*log(fc))
	coef_lambda<-1/coef_B
	########################################################################################
	# 04b. calc Saturated hydraulic conductivity Ksat [mm/hr] from calibrated Saxton (2006)
	######################################################################################## 
	ksat<-(6587*(sat-fc)^(3.347-coef_lambda ))
	
	#### Correction for peatlands very high SOM
	########################################################################################
	# 04. calc Saturated hydraulic conductivity Ksat [mm/hr] Balland et al. (2008)
	######################################################################################## 
	
	# L_10_Ksat<- -2.653985+3.092411*log10(dp-bd)+4.214627*sand
	# ksat_ball<-10^L_10_Ksat
	# # to mm/h
	# ksat_ball<-ksat_ball*10
	# 
	# 
	# if(!is.numeric(sand)){
	# 	ksat[ksat<=0] <- ksat_ball[ksat<=0]
	# }else{
	# 	ksat[!is.na(ksat) & ksat<=0]<-ksat_ball[!is.na(ksat) & ksat<=0]
	# }
	# 
	
		
	moist_fvol33init<-0.278*sand+0.034*clay+0.022*OM-0.018*(sand*OM)-0.027*(clay*OM)-0.584*(sand*clay)+0.078
	moist_fvol33<-moist_fvol33init+(0.636*moist_fvol33init-0.107)
	bub_init<--21.6*sand-27.93*clay-81.97*moist_fvol33+71.12*(sand*moist_fvol33)+8.29*(clay*moist_fvol33)+14.05*(sand*clay)+27.16
	bubbling_p<-bub_init+(0.02*bub_init^2-0.113*bub_init-0.7)
	# 101.97162129779 converts from KPa to mmH2O
	bubbling_p<-bubbling_p*-101.97162129779
	# If positive pressure: error in empirical fitting, water repelency or positive air entry pressure?? see  Wang et al., 2003, 10.1016/B978-0-444-51269-7.50009-6 :	#assume air entry pressure at the intercept of the log-log retention curve
	if(!is.numeric(sand)){
		bubbling_p[bubbling_p>0]<-coef_A[bubbling_p>0]*-101.97162129779
	}else{
		bubbling_p[!is.na(bubbling_p)& bubbling_p>0]<-coef_A[!is.na(bubbling_p)& bubbling_p>0]*-101.97162129779
	}
	
	#I still found around n=100 air entry pressures higher than 33kPa using the testing db (n=68567), assume 90% of saturation for those samples
	#bubpt[bubpt<=-3365.064]<-(33-(33*((sat-fc)/(0.9*sat-fc))))*101.97162129779
	########################################################################################
	# 05. calc shape parameters water retention for van Genuchten's curve Rawls (1985)
	######################################################################################## 
	sand<-sand*100
	clay<-clay*100
	silt<-100-sand-clay
	OM<-OM*100
		
	RES<--0.018+0.0009*sand+0.005*clay+0.029*sat -0.0002*clay^2-0.001*sand*sat-0.0002*clay^2*sat^2+0.0003*clay^2*sat -0.002*sat^2*clay
	RES[RES<0]<-0.0
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
	
	return(results)
}

julian_day <- function(y, m, i) {
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
	if (m <= 2) {
		y <- y - 1
		m <- m + 12
	}
	a <- floor(y/100)
	b <- 2 - a + floor(a/4)
	
	jde <- floor(365.25*(y + 4716)) + floor(30.6001*(m + 1)) + i + b - 1524.5
	return(jde)
}

dsin <- function(d) {
	# ************************************************************************
	# Name:     dsin
	# Inputs:   double (d), angle in degrees
	# Returns:  double, sine of angle
	# Features: This function calculates the sine of an angle (d) given
	#           in degrees.
	# Depends:  pir
	# ************************************************************************
	pir <- pi/180       # pi in radians
	sin(d*pir)
}

cum.interp<-function(x.months,y){
	# ************************************************************************
	# Name:     cum.interp, avg.interp
	# Inputs:   
	#               x.month ..... double, variable monthly value
	#               y ....... year
	# Returns:  double, daily values
	# Features: Interpolate monthly values to daily values
	# Depends:  julian day
	# ************************************************************************
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

frain_func<-function(tc,Tt,Tr,y,ny=NULL){
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
	
	# Features: calculates the snowfall occurrence probability
	# Ref:      Jennings, K.S., Winchell, T.S., Livneh, B., Molotch, N.P., 2018. Spatial variation of the 
	# 		  rain-snow temperature threshold across the Northern Hemisphere. Nat. Commun. 9, 1–9. doi:10.1038/s41467-018-03629-7
	# ************************************************************************
	# calibration set
	#p_snow<-1/(1+exp(-0.5827+1.319*as.numeric(tc)-as.numeric(elev)*4.18E-4-abs(as.numeric(lat))*1.140E-2))
	#all the data
	p_snow<-1/(1+exp(-0.4710405934+1.0473543991*as.numeric(tc)-as.numeric(elev)*0.0004596581-abs(as.numeric(lat))*0.0110592101))
	return(p_snow)
}

rspin_up <-function(lat,elev, sw_in, tc, pn, slop,asp, y,soil_data, Au,resolution) {
	
	# ************************************************************************
	# Name:     rspin_up
	# Inputs:   - vectors, forcing data
	#               lat ....................... double latitude (degrees)
	#               elev ...................... double, elevation (m)
	#               y ......................... int, year
	#               sw_in ..................... double, shortwave radiation
	#               tc ........................ double, mean daily air temperature (deg. C)
	#               pn ........................ double, precipitation (mm/d)
	#               soil_data ................. double, inputs for soil_hydro
	#               slope, asp, Au ............ double, terrain info
	# Returns:  list, daily fluxes and storages at steady state
	# Features: Wrapper of the c++ function, updates soil moisture and snow water equivalent until equilibrium
	# Depends:  soil_hydro, snowfall_prob, frain_func
	# ************************************************************************
	#### correct orientation slopes, from standard 0deg s north, in Allen, 2006 doi:10.1016/j.agrformet.2006.05.012 0deg is south!!!
	asp<-asp-180
	# get number of days in the year y
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	# interpolate monthly to daily
	if(length(sw_in)==12){sw_in<-avg.interp(sw_in,y)}
	if(length(tc)==12){tc<-avg.interp(tc,y)}
	if(length(pn)==12){pn<-cum.interp(pn,y)}
	# get soil hydrophysical characteristics
	soil_info<-soil_hydro(sand=soil_data[1],clay=soil_data[2],OM=soil_data[3],fgravel =soil_data[4] ,bd = soil_data[5])
	depth <- soil_data[6]
	SAT<-soil_info$SAT*depth*1000
	WP<-soil_info$WP*depth*1000
	FC<-soil_info$FC*depth*1000
	RES<-soil_info$RES*depth*1000
	lambda<-1/soil_info$B
	bub_press<-soil_info$bubbling_p
	if(length(Au)==1){
		#if there is no information on how many cell drain to this point, assume 3 sides of the octogonal cell
		ncellin<-3
		ncellout<-3
		soil_info<-c(SAT,WP,FC,soil_info$Ksat,lambda,depth,bub_press,RES,Au[1],resolution^2,ncellin,ncellout)
	}else{
		ncellin<-Au[2]
		ncellout<-Au[3]
		soil_info<-c(SAT,WP,FC,soil_info$Ksat,lambda,depth,bub_press,RES,Au[1],resolution^2,ncellin,ncellout)
	}
	

	# define snowfall occurrence:
	# 1. get snowfall probability of occurrence
	p_snow<-snowfall_prob(tc,lat,elev)
	# 2. get the threshold for snowfall occurrence
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

run_one_year <- function(lat,elev,slop,asp,sw_in, tc, pn, wn, y, snow,soil_data,Au,resolution,qin,td, snwage) {
	# ************************************************************************
	# Name:     run_one_year
	# Inputs:   - vectors, forcing data
	#               lat ....................... double latitude (degrees)
	#               elev ...................... double, elevation (m)
	#               y ......................... int, year
	#               sw_in ..................... double, shortwave radiation
	#               tc ........................ double, mean daily air temperature (deg. C)
	#               pn ........................ double, precipitation (mm/d)
	#               soil_data ................. double, inputs for soil_hydro
	#               slope, asp, Au ............ double, terrain info
	#               resolution ................ double, length of the cell (m) used to compute upslope area
	#               snow ...................... double, Snow water equivalent from the previous year
	#               wn ........................ double, soil water content from the previous year
	#               qin ....................... double, lateral flow input from upslope (mm/d)
	#               td ........................ double, time until drainage from upslope ceases (d)
	#               snwage .................... double, snow age (days)
	# Returns:  list, daily fluxes and storages for the year y
	# Features: Wrapper of the c++ function
	# Depends:  soil_hydro, snowfall_prob, frain_func
	# ************************************************************************
	#### correct orientation slopes, from standard 0deg s north, in Allen, 2006 doi:10.1016/j.agrformet.2006.05.012 0deg is south!!!
	asp<-asp-180
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	#interpolate if inputs are monthly
	if(length(sw_in)==12){sw_in<-avg.interp(sw_in,y)}
	if(length(tc)==12){tc<-avg.interp(tc,y)}
	if(length(pn)==12){pn<-cum.interp(pn,y)}
	if(length(wn)==12){wn<-avg.interp(wn,y)}
	#correct vector length for leap years
	if(length(wn)>length(pn)){
		wn<-wn[1:length(pn)]
		snow<-snow[1:length(pn)]
		qin<-qin[1:length(pn)]
		td<-td[1:length(pn)]
		snwage<-snwage[1:length(pn)]
	}else if(length(wn)<length(pn)){
		wn<-c(wn,wn[length(wn)])
		snow<-c(snow,snow[length(snow)])
		qin<-c(qin,qin[length(qin)])
		td<-c(td,td[length(td)])
		snwage<-c(snwage,snwage[length(snwage)])
	}
	
	# get soil hydrophysical characteristics
	soil_info<-soil_hydro(sand=soil_data[1],clay=soil_data[2],OM=soil_data[3],fgravel =soil_data[4] ,bd = soil_data[5])
	depth <- soil_data[6]
	SAT<-soil_info$SAT*depth*1000
	WP<-soil_info$WP*depth*1000
	FC<-soil_info$FC*depth*1000
	RES<-soil_info$RES*depth*1000
	lambda<-1/soil_info$B
	bub_press<-soil_info$bubbling_p
	if(length(Au)==1){
		#if there is no information on how many cell drain to this point, assume 3 sides of the octogonal cell
		ncellin<-3
		ncellout<-3
		soil_info<-c(SAT,WP,FC,soil_info$Ksat,lambda,depth,bub_press,RES,Au[1],resolution^2,ncellin,ncellout)
	}else{
		ncellin<-Au[2]
		ncellout<-Au[3]
		soil_info<-c(SAT,WP,FC,soil_info$Ksat,lambda,depth,bub_press,RES,Au[1],resolution^2,ncellin,ncellout)
	}
	# define snowfall occurrence:
	# 1. get snowfall probability of occurrence
	p_snow<-snowfall_prob(tc,lat,elev)
	# 3. get the threshold for snowfall occurrence
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
	# run splash for one year
	
	daily_totals<-my_splash$run_one_year(as.integer(ny), as.integer(y),as.numeric(sw_in),as.numeric(tc),as.numeric(pn),as.numeric(wn),slop,asp,as.numeric(snow),as.numeric(snowfall),soil_info,as.numeric(qin),as.numeric(td),as.numeric(snwage))
	return(daily_totals)
}

