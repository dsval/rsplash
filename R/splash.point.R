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
splash.point<-function(sw_in, tc, pn, lat,elev,slop=0,asp=0,soil_data,Au=0,resolution=250,time_index=NULL,monthly_out=FALSE,ts_out=TRUE,verbose=TRUE){
	###########################################################################
	# 010. Check time info
	###########################################################################	
	if(class(tc) %in% c('xts','zoo','ts')){
		time_index<-time(tc)
	}else if (!(class(tc)%in% c('xts','zoo','ts')) & !is.null(time_index)){
		if(verbose){message('Data is not a time-series object, getting time info from "time_index" ')}
		
		if(class(time_index)!="Date"){
			stop('time_index is not a "Date" class object')
		}
		
		
	}
	
	# else{
	# 	tc<-xts::xts(tc,time_index)
	# 	sw_in<-xts::xts(sw_in,time_index)
	# 	pn<-xts::xts(pn,time_index)
	# }
	###########################################################################
	# 01.Extract time info from the data
	###########################################################################
	# time frequency
	time.freq<-abs(as.numeric(time_index[1]-time_index[2], units = "days"))
	time_index_month<-seq(time_index[1],time_index[length(time_index)],'month')
	#### days of the year
	doys<-as.integer(format(time_index,'%j'))	
	# years
	yrs<-as.integer(format(time_index,'%Y'))
	y<-unique(yrs)
	# ndays in the year
	#ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	###########################################################################
	# 01.Get daily values if the inputs are monthly
	###########################################################################
	if(time.freq>20){
		time_index<-seq(as.Date(paste0(y[1],'-01-01')),as.Date(paste0(y[length(y)],'-12-31')),'day')
		time_index_month<-seq(time_index[1],time_index[length(time_index)],'month')
		#### days of the year
		doys<-as.integer(format(time_index,'%j'))	
		# years
		yrs<-as.integer(format(time_index,'%Y'))
		ndaypmonth<-as.data.frame(table(format(time_index,'%Y-%m')))
		####interpolations
		if (sum(!is.na(tc)) < 2) {
			tc<-rep(NA,length(time_index))} else{
			tc<-approx(time_index_month, tc, time_index, method = "linear", rule=2)$y
			}
		
		pn<-month2day_rain(pn,ndaypmonth$Freq)
		if (sum(!is.na(sw_in)) < 2) {
			sw_in<-rep(NA,length(time_index))} else{
			sw_in<-approx(time_index_month, sw_in, time_index, method = "linear", rule=2)$y
		}
		
	}
		
	
	###########################################################################
	# 03. initialize c++ program
	###########################################################################
	my_splash = new(SPLASH, lat, elev)
	###########################################################################
	# get soil hydrophysical characteristics
	###########################################################################
	soil_info<-soil_hydro(sand=soil_data[1],clay=soil_data[2],OM=soil_data[3],fgravel =soil_data[4] ,bd = soil_data[5])
	depth <- soil_data[6]
	SAT<-soil_info$SAT*depth*1000
	WP<-soil_info$WP*depth*1000
	FC<-soil_info$FC*depth*1000
	RES<-soil_info$RES*depth*1000
	Wmax<-soil_info$theta_c*depth*1000
	KWm<-soil_info$AWC*depth*1000
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
	###########################################################################
	# define snowfall occurrence:
	###########################################################################
	# 1. get snowfall probability of occurrence
	p_snow<-snowfall_prob(tc,lat,elev)
	# 2. get the threshold for snowfall occurrence
	Tt<-max(tc[p_snow>=0.5])
	# 3. get the fraction of precipitation falling as rain
	f_rain<-ifelse(p_snow>=0.5,frain_func(tc,Tt,13.3,time_index)[[1]],1)
	# define snowfall and rainfall:
	snowfall<-as.numeric(pn)*(1-f_rain)
	pn<-as.numeric(pn)*f_rain
	
	#### correct orientation slopes, from standard 0deg s north, in Allen, 2006 doi:10.1016/j.agrformet.2006.05.012 0deg is south!!!
	asp<-asp-180
	###########################################################################
	# 02. Get steady state
	###########################################################################
	# average daily for spin-up
	# sw_av<-tapply(sw_in,format(time(sw_in),"%j"),mean, na.rm=TRUE)
	# tc_av<-tapply(tc,format(time(sw_in),"%j"),mean, na.rm=TRUE)
	# pn_av<-tapply(pn,format(time(sw_in),"%j"),mean, na.rm=TRUE)
	
	# run spin up
	initial<-my_splash$spin_up(as.integer(365), as.integer(y[1]), as.numeric(sw_in[1:365]), as.numeric(tc[1:365]),as.numeric(pn[1:365]),slop,asp,as.numeric(snowfall[1:365]),soil_info)
	###########################################################################
	# run splash
	###########################################################################
	
	
	result<-my_splash$run_all(
		doys=doys,
		yrs=yrs,
		sw_in=as.numeric(sw_in),	# shortwave radiation W/m2
		tair=as.numeric(tc),		# air temperature C
		pn= pn,		# precipitation mm
		wn_last = initial$sm[length(initial$sm)],
		slop=slop,	# slope deg
		asp=asp,		# aspect deg
		snow_last=initial$snow[length(initial$snow)],
		snowfall = snowfall,
		soil_info=soil_info, 		# soil data: sand,clay,som in w/w %. Gravel v/v %, bulk density g/cm3, and depth to the bedrock (m)**
		qin_last=initial$qin[length(initial$qin)],
		td_last=initial$tdrain[length(initial$tdrain)],
		nds_last=initial$snwage[length(initial$snwage)])
	
		
	
	###########################################################################
	# 03. Start calculations if the  inputs are monthly 
	###########################################################################
	####################################################################################################
	####################################################################################################
	### subset the results
	result<-result[1:8]
	####################################################################################################
	# 4. compute soil moisture limitations
	####################################################################################################
	#get the wilting point and bucket size in mm
	#soil_water<-soil_hydro(sand=soil_data[1],clay=soil_data[2],OM=soil_data[3],fgravel =0.0 ,bd = soil_data[5])
	# wp<-(soil_water$WP/2)*soil_data[6]*1000
	# KWm<-(soil_water$FC-(soil_water$WP/2))*(soil_data[6]*1000)
	
	#get relative soil moisture limitation from 0.0 (at WP) to 1.0 (at FC)
	# soil_lim<-(result$wn-WP)/KWm
	# #adjust the boundaries, wn goes from ~WP to SAT
	# soil_lim[soil_lim<0]<-0.0
	# soil_lim[soil_lim>1]<-1.0
	# result$sm_lim<-soil_lim	
	soil_lim<-(result$wn-RES)/(Wmax-RES)
	#adjust the boundaries, wn goes from ~WP to SAT
	soil_lim[soil_lim<0]<-0.0
	soil_lim[soil_lim>1]<-1.0
	result$sm_lim<-soil_lim
	#result$sm_lim<-result$aet/result$pet	
	###########################################################################
	###########################################################################
	# 06. aggregate to monthly
	###########################################################################
	if(monthly_out){
		ind.mont<-format(time_index,'%Y-%m')
		###storages
		result[c('wn','snow','sm_lim')]<-lapply(result[c('wn','snow','sm_lim')],FUN=function(x){fastmatch::ctapply(x,INDEX=ind.mont,FUN=mean,na.rm=T)})
		result[c('ro','aet','pet','cond','bflow','netr')]<-lapply(result[c('ro','aet','pet','cond','bflow','netr')],FUN=function(x){fastmatch::ctapply(x,INDEX=ind.mont,FUN=sum,na.rm=T)})
		
				
	}
	
	if(ts_out & !monthly_out){
		result<-xts(base::as.data.frame(result),time_index)
	}else if (ts_out & monthly_out){
		result<-xts(base::as.data.frame(result),time_index_month)
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
	fsand<-sand/100
	fclay<-clay/100
	fOM<-OM/100
	fgravel<-fgravel/100
	########################################################################################
	# 02. calc bulk density [g/cm3] in case it is not provided, assumming 30 cm depth Balland et al. (2008)
	########################################################################################
	depth<-30
	dp<-1/((fOM/1.3)+((1-fOM)/2.65))
	
	if(!is.numeric(sand)){
		if(is.null(bd)){
			bd<-(1.5 + (dp-1.5-1.10*(1 - fclay))*(1-exp(-0.022*depth)))/(1+6.27*fOM)
		}
	}else{
		bd<-ifelse(is.na(bd),(1.5 + (dp-1.5-1.10*(1 - fclay))*(1-exp(-0.022*depth)))/(1+6.27*fOM),bd)
		
	}
	####error Fc calculations low bulk density. brute force low bd
	bd[bd<0.81]<-0.81
	########################################################################################
	# 03. calc volumetric water contents at saturation, 33kPa (fc) and 1500kPa (wp) Balland et al. (2008) - optimized
	######################################################################################## 
	# volumetric water content at saturation [m3/m3] first approx
	sat<-1-(bd/dp)
	# volumetric water content at 33kPa [m3/m3]
	fc<-(sat/bd)*(0.4760944 + (0.9402962 - 0.4760944)*fclay^0.5)*exp(-1*(0.05472678*fsand - 0.01* fOM)/(sat/bd))
	# volumetric water content at 1500kPa [m3/m3]
	wp_Ball<- fc*(0.2018522 + (0.7809203 - 0.2018522)*fclay^0.5) 
	########################################################################################
	# 04a. calc volumetric water contents at 1500kPa (wp) Sandoval(in prep.)
	######################################################################################## 
	# volumetric water content at 1500kPa [m3/m3]
	wp= -2.464e-05*sand+3.650e-03*clay+8.680e-03*OM+9.393e-03*bd  
	wp[!is.na(wp) & wp>=fc]<-wp_Ball[!is.na(wp) & wp>=fc]
	
	########################################################################################
	# 04b. calc volumetric water contents at saturation, 33kPa (fc) and 1500kPa (wp) Saxton and Rawls (2006)
	######################################################################################## 
	# FCinit<--0.251*sand+0.195*clay+0.011*OM+0.006*(sand*OM)-0.027*(clay*OM)+0.452*(sand*clay)+0.299
	# FC_fvol<-FCinit+(1.283*FCinit^2-0.374*FCinit-0.015)
	# fc[fc<0 | fc>=sat]<-FC_fvol[fc<0 | fc>=sat]
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
	#wp<- fc*(0.2018522 + (0.7809203 - 0.2018522)*clay^0.5) 
	
	########################################################################################
	# 05. calc shape parameters water retention Brooks and Corey curve Saxton and Rawls (2006)
	######################################################################################## 
	# get parameters for BC eqn form Saxton 2006
	coef_B<-(log(1500)-log(33))/(log(fc)-log(wp))
	coef_A<-exp(log(33)+coef_B*log(fc))
	coef_lambda<-1/coef_B
	########################################################################################
	# 06a. calc theta crit aasuming z 2m
	######################################################################################## 
	#coeff_c = 1000.0/(pw*G);
	coeff_c = 1000.0/(997*9.80665);
	theta_c<-(coeff_c*coef_A/2.0)^(1/(1+coef_B))
	########################################################################################
	# 06a. calc theta crit aasuming z 1m
	######################################################################################## 
	#coeff_c = 1000.0/(pw*G);
	#coeff_c = 1000.0/(997*9.80665);
	#theta_c<-(coeff_c*coef_A/1.0)^(1/(1+coef_B))
	########################################################################################
	# 07. residual water content lm Sandoval(in prep.)
	######################################################################################## 
	#theta_r<-0.4068341*wp+0.0050417*coef_lambda+0.0143244*bd
	theta_r<-0.53097*wp
	########################################################################################
	# 08a. calc Saturated hydraulic conductivity Ksat [mm/hr] from calibrated Saxton and Rawls (2006)
	######################################################################################## 
	#ksat<-(4743.514223*(sat-fc)^(2.737082+0.240194*coef_lambda ))
	#ksat<-(6587*(sat-fc)^(3.347-coef_lambda ))
	

########################################################################################
# 09. calc Air entry pressure [mm] from calibrated Saxton and Rawls (2006)
######################################################################################## 

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
# 08b. calc Saturated hydraulic conductivity Ksat [mm/hr] from rsplash. Sandoval(in prep.)
######################################################################################## 

sat<-sat*(1-fgravel)
fc<-fc*(1-fgravel)
wp<-wp*(1-fgravel)
#ksat <- 1200/(1 + exp(7.5623 -17.1442 * (sat-fc) -0.9141 * bd + 0.7076 * coef_lambda))
ksat <- 4.726e+02/(1 + exp(2.670*bd+92.1*theta_r+12.38*(fc-wp)-9.409 * 
		theta_c -1.416e+01 * (sat-fc) + 6.499e-15 * bubbling_p ))
#coeff_ksat<-c(ksmax=857.48454,k2=-2.70927 ,k3=3.62264,k4=7.33398,k5=-8.11795, k6=18.75552, k7=1.03319 )

# 	ksat<-coeff_ksat['ksmax']/(1 + exp(coeff_ksat['k2'] * fsand + 
# 		coeff_ksat['k3'] * bd + 
# 		coeff_ksat['k4'] * fclay + 
# 		coeff_ksat['k5'] *	(sat-fc) + 
# 		coeff_ksat['k6'] * fOM + 
# 		coeff_ksat['k7'] * coef_lambda)
# )

########################################################################################
# 05. calc shape parameters water retention for van Genuchten's curve Rawls (1985)
######################################################################################## 
silt<-100-sand-clay

# RES<--0.018+0.0009*sand+0.005*clay+0.029*sat -0.0002*clay^2-0.001*sand*sat-0.0002*clay^2*sat^2+0.0003*clay^2*sat -0.002*sat^2*clay
# RES[RES<0]<-0.0
# parameters for van Genutchen eqn
topsoil<-1

alpha<-exp(-14.96 + 0.03135*clay + 0.0351*silt + 0.646*OM +15.29*dp - 0.192*topsoil -4.671*dp^2- 0.000781*clay^2 - 0.00687*OM^2 + 0.0449/OM + 0.0663*log(silt) + 0.1482*log(OM) - 0.04546*dp *silt - 0.4852*dp*OM + 0.00673*topsoil*clay)

n<-1.0+exp(-25.23 - 0.02195*clay + 0.0074*silt - 0.1940*OM + 45.5*dp - 7.24*dp^2 +0.0003658*clay^2 + 0.002885*OM^2 -12.81/dp - 0.1524/silt - 0.01958/OM - 0.2876*log(silt) - 0.0709*log(OM) -44.6*log(dp) - 0.02264*dp*clay + 0.0896*dp*OM +0.00718*topsoil*clay)

m<-1-(1/n)


results$SAT<-sat
results$FC<-fc
results$WP<-wp
results$bd<-bd
results$AWC<-(fc-wp)
results$Ksat<-ksat
results$A<-coef_A
results$B<-coef_B
results$theta_c<-theta_c
results$RES<-theta_r*(1-fgravel)
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



month2day_rain<-function(pn_monthly,ndaypmonth){
	# ************************************************************************
	# Name:     month2day_rain
	# Inputs:   
	#               pn_monthly ..... double, variable monthly value
	#               y ....... year
	# Returns:  double, daily values
	# Features: Generate daily synthetic values of precipitation
	# Depends:  base R
	# Ref.  Geng, et al., (1986) doi:10.1016/0168-1923(86)90014-6
	# ************************************************************************
	###get time info
	
	### get parameters for a Gamma distribution as in Geng, et al., (1985, 1986)
	geo_mean<-exp(mean(log(pn_monthly[pn_monthly>0.0]),na.rm=T))
	mean_mean<-mean(pn_monthly,na.rm=T)
	if(!is.na(mean_mean) & mean_mean>0){
		Y<-log(mean_mean/geo_mean)
		
		if(Y>0.0 & Y < 0.5772){
			alph<-(0.5000876 +0.16488552*Y-0.0544274*Y^2)/Y
		}else if (Y>=0.5772){
			alph<-(8.898919+9.059950*Y+0.9775373*Y^2)/(Y*(17.79728+11.968477*Y+Y^2))
		}else{
			alph<-1
		}
		
		bet<-mean_mean/alph
	}else{
		alph<-NA;bet<-NA
	}
	
	
	
	#####################################################################################
	
	get_pn_day<-function(ndays,pmonth,alph,bet){
		if(is.na(pmonth)){
			return(rep(NA,ndays))
		}else if(pmonth==0){
			return(rep(0.0,ndays))
		}else{
			##daily precipitatio uncorrected
			pday<-rgamma(ndays,alph,bet)
			##correction factor
			fac=pmonth/sum(pday)
			##daily corrected precipitation
			return(pday*fac)
		}
		
	}
	
	
	daily_pn<-mapply(get_pn_day,ndaypmonth,pn_monthly,MoreArgs = list(alph,bet),SIMPLIFY = F)
	
	do.call(c,daily_pn)
}




frain_func<-function(tc,Tt,Tr,time_index){
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
	
	# if(length(tc)>1){
	# 	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	# 	ndaysmonth<-rep(NA,12)
	# 	for(i in 1: 12){ndaysmonth[i]<-julian_day(y,i+1,1)-julian_day(y,i,1)}
	# 	m_ind<-rep(1,ndaysmonth[1])
	# 	for(i in 2:12){m_ind<-c(m_ind,rep(i,ndaysmonth[i]))}	
	# }else{
	# 	m_ind<-as.Date(paste(y,ny,sep="-"),format="%Y-%j")
	# 	m_ind<-format(m_ind,"%m")
	# 	m_ind<-as.numeric(m_ind)
	# }	
	#m_ind<-as.Date(paste(y,ny,sep="-"),format="%Y-%j")
	m_ind<-format(time_index,"%m")
	m_ind<-as.numeric(m_ind)
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

