#' unsSWC.grid
#'
#' Estimates water content in the unsaturated part of the profile using an analytical integral of the eqn. described by Brooks and Corey 1964 and calculates depth to water table if there is any (*experimental)
#' @param  soil_data, Raster* object with the layers organized as sand(perc),clay(perc),organic matter(perc),coarse-fragments-fraction(perc), bulk density(g cm-3)
#' @param  uns_depth, maximum depth to calculate the accumulated water content in (m)
#' @param  wn, Raster* object, total water content of the profile from surface to bedrock (m)
#' @param  outdir, directory where to save the files
#' @import raster
#' @keywords splashTools
#' @export
#' @examples
#' unSWC()

unSWC.grid<-function(soil_data,uns_depth,wn,d_name ,outdir=getwd()){
#########################################################################
######################### 1.0 define functions ##############################
#########################################################################
	ztime<-getZ(wn)
	time.freq<-ztime[2]-ztime[1]
		
	if(abs(as.numeric(time.freq, units = "days"))<2){
		time_text<-'daily'
		}else{
		time_text<-'monthly'
	}
	
	
	
	

	UnsWater<-function(psi_m,uns_depth,theta_r,theta_s,bub_press,lambda,depth){
		#-----------------------------------------------------------------------
		# Input:    - float, Matric potential (Psi_m), mm
		#           - float, Depth to where to calculate water content (uns_depth), m
		#           - float, residual moisture (theta_r), m3 m-3
		#           - float, saturated moisture (theha_s), m3 m-3
		#           - float, bubbling/air entry pressure (bub_press), m3 m-3
		#           - float, slope of the log curve thetha vs Psi_m (lambda), m3 m-3
		# Output:   float, integrated water from surface to depth z_uns, mm
		# Features: Estimates water content using an analytca integral of the profile described by Brooks and Corey 1964
		#          
		# Ref:      Brooks, R.H., Corey, A.T., 1964. Hydraulic properties of porous media. 
		#   		  Hydrology Papers No 17. Colorado State University. doi:10.13031/2013.40684
		#           
		#-----------------------------------------------------------------------
		#1.Calc wtd
		wtdini<-(bub_press-psi_m)/1000
		wtd<-ifelse(wtdini>depth,depth,ifelse(wtdini<0,0,wtdini))
		#2. calc unsaturated depth
		z_uns<-ifelse(wtd<=uns_depth,wtd*1000,uns_depth*1000)
		# z_uns<-uns_depth*1000
		# 3. calc analytical integral from 0 to z_uns
		w_uns_z<- theta_r*z_uns+(((psi_m+z_uns)*(theta_r-theta_s)*(bub_press/(psi_m+z_uns))^lambda)/(lambda-1))
		w_uns_0<- theta_r*0+(((psi_m+0)*(theta_r-theta_s)*(bub_press/(psi_m+0))^lambda)/(lambda-1))
		w_uns<-w_uns_z-w_uns_0
		#4. 
		sat_swc<-ifelse(wtd<=uns_depth,theta_s*(uns_depth-wtd)*1000,0)
		
		swc_tot<-w_uns+sat_swc
		# w_uns[psi_m>=bub_press]<-theta_s*z_uns
		return(swc_tot)
		# return(list(wn_z=swc_tot,wtd=wtd))
		
	}
	
	# wn<-simdf[,1]
	#########################################################################
	######################### 2. get soil hydrophysics ##############################
	soil_info<-soil_hydro(sand=soil_data[[1]],clay=soil_data[[2]],OM=soil_data[[3]],fgravel =soil_data[[4]]*0 ,bd = soil_data[[5]])
	## assume 0.0 as residula water 
	theta_s<-soil_info$SAT
	theta_r<-soil_info$RES
	z_uns <- uns_depth
	lambda<-1/soil_info$B
	bub_press<-soil_info$bubbling_p
	#theta_i<-(wn/(soil_data[[6]]*1000))#*(1-soil_data[[4]]/100)
	#########################################################################
	######################### 2. calc mean theta ##############################
	setwd(outdir)
	### get time info
	y=as.integer(unique(format(getZ(wn),'%Y')))
	calc_thetai<-function(w,s,r,d){
		theta_o<-(w/(d*1000))
		w<-ifelse(theta_o>=s,s-0.0001,ifelse(theta_o<=r,r+0.0001,theta_o))
		
		# i[i>=s]<-s-0.05
		# i[i<=r]<-r+0.01
		w
	}
	theta_i<-overlay(wn,theta_s,theta_r,soil_data[[6]],fun=calc_thetai,filename=paste0(outdir,"/",'SPLASH_v2.0_',"theta_mean",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="theta", varunit="m3/m3",longname="volumetric soil moisture", xname="lon", yname="lat", zname="time")
	
	# writeRaster(result.all$wn,paste0(outdir,"/",y[1],"_",y[length(y)],".","theta_top",".","nc"),format="CDF",overwrite=TRUE,varname="theta", varunit="m3/m3",longname="volumetric soil moisture", xname="lon", yname="lat", zname="time")
	#########################################################################
	######################### 4. calc mean matric potential ##############################
	psi_m = bub_press/((((theta_i-theta_r)/(theta_s-theta_r)))^(1/lambda));
	#########################################################################
	######################### 5. calculate the water table depth (m)s ##############################
	# 
	calcwtd<-function(psi_m,totdepth,bub_press){
		wtdini<-(bub_press-psi_m)/1000
		wtd<-ifelse(wtdini>totdepth,totdepth,ifelse(wtdini<0,0,wtdini))
		# # wtd[wtd>totdepth]<-totdepth
		# # wtd[wtd<0]<-0
		wtd
		# wtdini
	}
	
	wtd<-overlay(psi_m,soil_data[[6]],bub_press,fun=calcwtd,filename=paste0(outdir,"/",'SPLASH_v2.0_',"wtd",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="wtd", varunit="m",longname="water table depth", xname="lon", yname="lat", zname="time")
	
	# # Find where the water table is whitin the root zone
	# shallow_wtd<-wtd<=uns_depth
	# # calculate the swc of the saturated part whitin the zoot zone
	# sat_swc<-(theta_s*(uns_depth-wtd)*1000)*shallow_wtd
	#########################################################################
	######################### 6. update unsaturated zone ##############################
	# update unsaturated zone
	w_z<-overlay(psi_m,z_uns,theta_r,theta_s,bub_press,lambda,uns_depth,fun=UnsWater,filename=paste0(outdir,"/",'SPLASH_v2.0_',"swc_top_",d_name,'_m','_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="swc", varunit="mm",longname="soil water content", xname="lon", yname="lat", zname="time")
	#########################################################################
	######################### 7. get effective saturation at the top ##############################
	#soil_info<-soil_hydro(sand=soil_data[[1]],clay=soil_data[[2]],OM=soil_data[[3]],fgravel =soil_data[[4]],bd = soil_data[[5]])
	calc_Se=function(w_z,uns_depth,theta_s){
		theta_i_top=w_z/(uns_depth*1000)
		Se=theta_i_top/theta_s
		Se=ifelse(Se>1,1,ifelse(Se<0,0,Se))
		Se
	}
	Se=overlay(w_z,uns_depth,theta_s,fun=calc_Se,filename=paste0(outdir,"/",'SPLASH_v2.0_',"Se_top_",d_name,'_m','_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="Se", varunit="fraction",longname="effective saturation or water filled porosity", xname="lon", yname="lat", zname="time")
	
	
	# psi_m<-as.numeric(psi_m)
	# uns_wn<-UnsWater(psi_m,z_uns,theta_r,theta_s,bub_press,lambda)
	# # uns_wn<-overlay(theta_i,theta_s,theta_r,uns_wn,z_uns,fun=cor_t_i)
	# #uns_theta<-uns_wn/(uns_depth*1000)
	# # uns_theta<-xts(uns_theta,time(wn))
	# # uns_wn<-xts(uns_wn,time(wn))
	# wtd<-(bub_press-psi_m)/1000
	# wtd[wtd>soil_data[[6]]]<-soil_data[[6]][wtd>soil_data[[6]]]
	return(list(w_z=w_z,wtd=wtd,Se=Se))
	# merge.xts(uns_wn,wtd)
	
}
