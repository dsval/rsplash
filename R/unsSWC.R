#' unsSWC
#'
#' Estimates water content using an analytca integral of the profile described by Brooks and Corey 1964 and calculates depth to water table if there is any
#' @param  soil_data, vector as retrieved from splashTools::getSoilNRCS() or splashTools::getSoilSite()
#' @param  uns_depth, maximum depth to calculate the accumulated water content in (m)
#' @param  wn, xts object, total water content of the profile from surface to bedrock (m)
#' @param  ouputdir, directory where to save the files
#' @import xts
#' @keywords splashTools
#' @export
#' @examples
#' unSWC()

unSWC<-function(soil_data,uns_depth,wn){
	# test
	# soil_data=soil_sm[[1]];uns_depth=sites_sm[1]*-1;wn=sm_sim$`AR-SLu`[,1]
	# end testing
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
		
		swc_tot<-(w_uns+sat_swc)
		# w_uns[psi_m>=bub_press]<-theta_s*z_uns
		# return(swc_tot)
		return(list(wn_z=swc_tot,wtd=wtd))
		
	}
	### get volumetric water content without stoniness
	soil_info<-soil_hydro(sand=soil_data[1],clay=soil_data[2],OM=soil_data[3],fgravel =0.0 ,bd = soil_data[5])
	wn<-wn*(100+soil_data[4])/100
	
	#########################################################################################
	theta_s<-as.numeric(soil_info$SAT)
	theta_r<-as.numeric(soil_info$RES)
	#theta_r<-0.0
	lambda<-as.numeric(1/soil_info$B)
	bub_press<-as.numeric(soil_info$bubbling_p)
	theta_i<-wn/(soil_data[6]*1000)
	theta_i[theta_i>=theta_s]<-theta_s-0.0001
	theta_i[theta_i<=theta_r]<-theta_r+0.0001
	psi_m = bub_press/((((theta_i-theta_r)/(theta_s-theta_r)))^(1/lambda));
	psi_m<-as.numeric(psi_m)
	#use Saxton and Rawls 2006
	#psi_m<-soil_info$A*(theta_i)^(-1*soil_info$B)
	#psi_m<-as.numeric(psi_m)*-101.97162129779
	uns_wn<-UnsWater(psi_m,uns_depth,theta_r,theta_s,bub_press,lambda,soil_data[6])
	###error at very low swc, super high potential, swc getsclose to 0
	uns_wn$wn_z[uns_wn$wn_z<=0]<-theta_r*(uns_depth*1000)
	
	if(is.xts(wn)){
		uns_wn<-mapply(FUN=xts,uns_wn,MoreArgs = list(order.by=time(wn)),SIMPLIFY = F)
	}
	
	uns_theta<-(uns_wn[[1]]/(uns_depth*1000))
	#### correct swc for stoniness at the top
	uns_wn[[1]]<-uns_wn[[1]]*(1-(soil_data[4]/100))
		
	result<-uns_wn
	result$theta_i<-uns_theta
		
	
	result
	
}
