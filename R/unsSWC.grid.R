#' unsSWC.grid
#'
#' Estimates water content using an analytca integral of the profile described by Brooks and Corey 1964 and calculates depth to water table if there is any
#' @param  soil_data, Raster* object with the layers organized as sand(perc),clay(perc),organic matter(perc),coarse-fragments-fraction(perc), bulk density(g cm-3)
#' @param  uns_depth, maximum depth to calculate the accumulated water content in (m)
#' @param  wn, Raster* object, total water content of the profile from surface to bedrock (m)
#' @param  ouputdir, directory where to save the files
#' @import raster
#' @keywords splashTools
#' @export
#' @examples
#' unSWC()

unSWC.grid<-function(soil_data,uns_depth,wn){
	# test
	# soil_data=soil_snow[[1]];uns_depth=depth_sm_snow[[1]];wn=sims$`SNTL:1243`[,1]
	# soil_data=soil50km;uns_depth=rootd;wn=wnd
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
		
		swc_tot<-w_uns+sat_swc
		# w_uns[psi_m>=bub_press]<-theta_s*z_uns
		return(swc_tot)
		# return(list(wn_z=swc_tot,wtd=wtd))
		
	}
	
	# wn<-simdf[,1]
	
	soil_info<-soil_hydro(sand=soil_data[[1]],clay=soil_data[[2]],OM=soil_data[[3]],fgravel =soil_data[[4]] ,bd = soil_data[[5]])
	# SAT<-soil_info$SAT*(1-soil_data[[4]]/100)*soil_data[[6]]*1000
	# RES<-soil_info$WP*(1-soil_data[[4]]/100)*soil_data[[6]]*1000
	# theta_s<-SAT/(soil_data[[6]]*1000)
	# theta_r<-RES/(soil_data[[6]]*1000)
	theta_s<-soil_info$SAT
	theta_r<-soil_info$WP
	z_uns <- uns_depth
	lambda<-1/soil_info$B
	bub_press<-soil_info$bubbling_p
	#theta_i<-(wn/(soil_data[[6]]*1000))#*(1-soil_data[[4]]/100)
	calc_thetai<-function(w,s,r,d){
		theta_o<-(w/(d*1000))
		w<-ifelse(theta_o>=s,s-0.01,ifelse(theta_o<=r,r+0.001,theta_o))
		
		# i[i>=s]<-s-0.05
		# i[i<=r]<-r+0.01
		w
	}
	theta_i<-overlay(wn,theta_s,theta_r,soil_data[[6]],fun=calc_thetai)
	
	psi_m = bub_press/((((theta_i-theta_r)/(theta_s-theta_r)))^(1/lambda));
	# 
	# calculate the water table depth (m)
	calcwtd<-function(psi_m,totdepth,bub_press){
		wtdini<-(bub_press-psi_m)/1000
		wtd<-ifelse(wtdini>totdepth,totdepth,ifelse(wtdini<0,0,wtdini))
		# # wtd[wtd>totdepth]<-totdepth
		# # wtd[wtd<0]<-0
		wtd
		# wtdini
	}
	
	wtd<-overlay(psi_m,soil_data[[6]],bub_press,fun=calcwtd)
	# # Find where the water table is whitin the root zone
	# shallow_wtd<-wtd<=uns_depth
	# # calculate the swc of the saturated part whitin the zoot zone
	# sat_swc<-(theta_s*(uns_depth-wtd)*1000)*shallow_wtd
	# update unsaturated zone
	w_z<-overlay(psi_m,z_uns,theta_r,theta_s,bub_press,lambda,uns_depth,fun=UnsWater)
	
	
	# psi_m<-as.numeric(psi_m)
	# uns_wn<-UnsWater(psi_m,z_uns,theta_r,theta_s,bub_press,lambda)
	# # uns_wn<-overlay(theta_i,theta_s,theta_r,uns_wn,z_uns,fun=cor_t_i)
	# #uns_theta<-uns_wn/(uns_depth*1000)
	# # uns_theta<-xts(uns_theta,time(wn))
	# # uns_wn<-xts(uns_wn,time(wn))
	# wtd<-(bub_press-psi_m)/1000
	# wtd[wtd>soil_data[[6]]]<-soil_data[[6]][wtd>soil_data[[6]]]
	return(list(w_z=w_z,wtd=wtd))
	# merge.xts(uns_wn,wtd)
	
}
