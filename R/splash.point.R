#' splash.point
#'
#' Apply splash algorithm
#' @param   sw_in, lon
#' @param   tc, lon
#' @param   pn, lon
#' @param   elev, lon
#' @return a matrix xts type
#' @import Rcpp 
#' @import xts
#' @keywords splash
#' @export
#' @examples
#' splash.grid()
splash.point<-function(sw_in, tc, pn, lat,elev,slop,asp,soil_data,Au,resolution){
	
	# require(xts)
	# Extract time info from data
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
			
			result<-run_one_year(lat,elev,slop,asp,sw_in, tc, pn,initial$sm, y[1], initial$snow,soil_data,Au,resolution,initial$qin,initial$tdrain)
			# result<-xts(result,ztime)
			result<-do.call(cbind,result)
		}
		else if(length(y)>1){
			
			end<-cumsum(ny)
			start<-end+1
			result<-list()
			sw_av<-tapply(sw_in,format(time(sw_in),"%j"),mean, na.rm=TRUE)
			tc_av<-tapply(tc,format(time(sw_in),"%j"),mean, na.rm=TRUE)
			pn_av<-tapply(pn,format(time(sw_in),"%j"),mean, na.rm=TRUE)
			# initial<-rspin_up(lat,elev, sw_in[1:ny[1]], tc[1:ny[1]], pn[1:ny[1]], slop,asp, y[1],soil_data,Au,resolution)
			initial<-rspin_up(lat,elev, sw_av, tc_av, pn_av, slop,asp, y[1],soil_data,Au,resolution)
			result[[1]]<-run_one_year(lat,elev,slop,asp,sw_in[1:ny[1]], tc[1:ny[1]],  pn[1:ny[1]],initial$sm, y[1], initial$snow,
				soil_data,Au,resolution,initial$qin,initial$tdrain)
			
			for (i in 2:length(y)){
				
				stidx<-i-1
				# correct for leap years	
				result[[i]]<-run_one_year(lat,elev,slop,asp, sw_in[start[stidx]:end[i]], tc[start[stidx]:end[i]], pn[start[stidx]:end[i]],
					result[[stidx]]$wn,y[i],result[[stidx]]$snow,soil_data,Au,resolution,result[[stidx]]$qin_prev,result[[stidx]]$tdrain)
			}
			result<-lapply(result,FUN=base::as.data.frame)
			result<-do.call(rbind,result)
			
		}
		result<-xts(result,ztime)
	}
	###########################################################################
	# 03. Start the calculations for Monthly inputs
	###########################################################################
	if (time.freq>20){		
		ztime.days<-seq(as.Date(paste(y[1],1,sep="-"),format="%Y-%j"),as.Date(paste(y[length(y)],ny[length(y)],sep="-"),format="%Y-%j"), by="day")
		if (length(y)==1){
			initial<-rspin_up(lat,elev, sw_in, tc, pn, slop,asp, y[1],soil_data,Au,resolution)
			
			result<-run_one_year(lat,elev,slop,asp,sw_in, tc, pn,initial$sm, y[1], initial$snow,soil_data,Au,resolution,initial$qin,initial$tdrain)
			# result<-xts(result,ztime)
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
			# initial<-rspin_up(lat,elev, sw_in[1:ny[1]], tc[1:ny[1]], pn[1:ny[1]], slop,asp, y[1],soil_data,Au,resolution)
			initial<-rspin_up(lat,elev, sw_av, tc_av, pn_av, slop,asp, y[1],soil_data,Au,resolution)
			result[[1]]<-run_one_year(lat,elev,slop,asp,sw_in[1:end[1]], tc[1:end[1]],  pn[1:end[1]],initial$sm, y[1], initial$snow,
				soil_data,Au,resolution,initial$qin,initial$tdrain)
			
			for (i in 2:length(y)){
				
				stidx<-i-1
				# correct for leap years	
				result[[i]]<-run_one_year(lat,elev,slop,asp, sw_in[start[i]:end[i]], tc[start[i]:end[i]], pn[start[i]:end[i]],
					result[[stidx]]$wn,y[i],result[[stidx]]$snow,soil_data,Au,resolution,result[[stidx]]$qin_prev,result[[stidx]]$tdrain)
			}
			result<-lapply(result,FUN=base::as.data.frame)
			result<-do.call(rbind,result)
			
		}
		result<-xts(result,ztime.days)
	}
	
			
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
	
	if(!is.numeric(sand)){
		if(is.null(bd)){
			bd<-(1.5 + (dp-1.5-1.10*(1 - clay))*(1-exp(-0.022*depth)))/(1+6.27*OM)
		} 
	}else{
		if(is.na(bd)){
			bd<-(1.5 + (dp-1.5-1.10*(1 - clay))*(1-exp(-0.022*depth)))/(1+6.27*OM)
		} 
	}
	
	
	
	
	sat<-1-(bd/dp)
	
	# fc<-(sat/bd)*(0.565 + (0.991 - 0.565)*clay^0.5)*exp(-(0.103*sand - 0.785* OM)/(sat/bd))
	fc<-(sat/bd)*(0.3366685 + (1.417544 - 0.3366685)*clay^0.5)*exp(-1*(0.03320495*sand - 0.2755312* OM)/(sat/bd))
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
	
	# results$VG_m<-m
	
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
	
	# Features: calculate the probability of snowfall occurrence probability, if >0.5, snowfall will occur
	# Ref:      Jennings, K.S., Winchell, T.S., Livneh, B., Molotch, N.P., 2018. Spatial variation of the 
	# 		  rain-snow temperature threshold across the Northern Hemisphere. Nat. Commun. 9, 1–9. doi:10.1038/s41467-018-03629-7
	# ************************************************************************
	p_snow<-1/(1+exp(-0.5827+1.319*as.numeric(tc)-as.numeric(elev)*4.18E-4-abs(as.numeric(lat))*1.140E-2))
	return(p_snow)
}

rspin_up <-function(lat,elev, sw_in, tc, pn, slop,asp, y,soil_data, Au,resolution) {
	
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
		ncellin<-2
		ncellout<-1
	}else{
		ncellin<-Au[2]
		ncellout<-Au[3]
	}
	
	soil_info<-c(SAT,WP,FC,soil_info$Ksat,lambda,depth,bub_press,RES,Au[1],resolution^2,ncellin,ncellout)
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
	# result<-my_splash$spin_up(as.integer(ny), as.integer(y), as.numeric(sw_av), as.numeric(tc_av),as.numeric(pn_av),slop,asp,as.numeric(snowfall),soil_info)
	
	return(result)
}

run_one_year <- function(lat,elev,slop,asp,sw_in, tc, pn, wn, y, snow,soil_data,Au,resolution,qin,td) {
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
	SAT<-soil_info$SAT*depth*1000
	WP<-soil_info$WP*depth*1000
	FC<-soil_info$FC*depth*1000
	RES<-soil_info$RES*depth*1000
	lambda<-1/soil_info$B
	bub_press<-soil_info$bubbling_p
	if(length(Au)==1){
		ncellin<-2
		ncellout<-1
	}else{
		ncellin<-Au[2]
		ncellout<-Au[3]
	}
	soil_info<-c(SAT,WP,FC,soil_info$Ksat,lambda,depth,bub_press,RES,Au[1],resolution^2,ncellin,ncellout)
	
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
	
	daily_totals<-my_splash$run_one_year(as.integer(ny), as.integer(y),as.numeric(sw_in),as.numeric(tc),as.numeric(pn),as.numeric(wn),slop,asp,as.numeric(snow),as.numeric(snowfall),soil_info,as.numeric(qin),as.numeric(td))
	return(daily_totals)
}

