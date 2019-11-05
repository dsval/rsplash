#' upslope_area
#'
#' Computes the upslope area in m2 using topmodel for small areas or Taudem for large
#' @param   sw_in, lon
#' @param   tc, lon
#' @param   pn, lon
#' @param   elev, lon
#' @return a matrix xts type
#' @import topmodel
#' @import raster 
#' @keywords splash
#' @export
#' @examples
#' splash.grid()
upslope_area<-function(dem){
	# require(raster)
	# setwd(tmpdir)
	# require(topmodel)
	# rasterOptions(maxmemory=3e7, timer=FALSE, tmptime = 24, chunksize = 3e7,todisk=FALSE, overwrite=TRUE)
	resolution<-sqrt(cellStats(area(dem), stat='mean', na.rm=TRUE))*1000
	elev<-sinkfill(raster::as.matrix(dem),resolution,1)
	elev<-raster(elev)
	extent(elev)<-extent(dem)
	crs(elev)<-crs(dem)
	
	# 01. Calculate spatial distributed variables
	# calculate watershed area
	areacatch<-topidx(raster::as.matrix(elev), resolution= resolution)$area
	areacatch<-raster(areacatch)
	extent(areacatch)<-extent(elev)
	crs(areacatch)<-crs(elev)
	areacatch<-writeRaster(areacatch,"areacatch.grd",overwrite=TRUE)
	flowdir<-terrain(elev,opt='flowdir')
	ncellin<-ncellflow(flowdir,inout='in',met='top',filename="ncellin.grd",overwrite=TRUE)
	ncellout<-ncellflow(flowdir,inout='out',met='top',filename="ncellout.grd",overwrite=TRUE)
	rm(elev)
	gc()
	# return(areacatch)
	return(stack(areacatch,ncellin,ncellout))
	
}
upslope_areav2<-function(dem,...){
	# returns upslope area in square meters
	# require(raster)
	# Set working directory to your location
	# setwd(tmpdir)
	writeRaster(dem,"rawdem.tif",format="GTiff", overwrite=TRUE)
	
	if(... == 'MPI'){
		# Pitremove
		system("pitremove -z rawdem.tif -fel dem_nopit.tif",show.output.on.console=F,invisible=F)
		
		# D8 flow directions
		system("d8flowdir -p dem_p.tif -sd8 dem_sd8.tif -fel dem_nopit.tif",show.output.on.console=F,invisible=F)
		# Contributing area
		system("aread8 -p dem_p.tif -ad8 dem_a_ac.tif -nc",show.output.on.console=F,invisible=F)
	}else{
		# Pitremove
		system("mpiexec pitremove -z rawdem.tif -fel dem_nopit.tif",show.output.on.console=F,invisible=F)
		
		# D8 flow directions
		system("mpiexec d8flowdir -p dem_p.tif -sd8 dem_sd8.tif -fel dem_nopit.tif",show.output.on.console=F,invisible=F)
		# Contributing area
		system("mpiexec aread8 -p dem_p.tif -ad8 dem_a_ac.tif -nc",show.output.on.console=F,invisible=F)
	}
	
	
	ups_ncell<-raster("dem_a_ac.tif")
	flowdir<-raster("dem_p.tif")
	area_p_cell<-area(ups_ncell)
	ups_area<-overlay(ups_ncell, area_p_cell, fun=function(x,y){(x*y*1e6)},filename="areacatch.grd",overwrite=TRUE)
	ncellin<-ncellflow(flowdir,inout='in',met='tau',filename="ncellin.grd",overwrite=TRUE)
	ncellout<-ncellflow(flowdir,inout='out',met='tau',filename="ncellout.grd",overwrite=TRUE)
	return(stack(ups_area,ncellin,ncellout))
	# return(ups_area)
	
}


ncellflow<-function(flowdir,inout='in',met='top', ...){
	if(inout=='in'&& met=='top'){
		flow<-matrix(c(2,1,128,4,0,64,8,16,32), nrow=3)
	}else if (inout=='out'&& met=='top'){
		flow<-matrix(rev(c(2,1,128,4,0,64,8,16,32)), nrow=3)
	}else if (inout=='in'&& met=='tau'){
		flow<-matrix(c(8,1,2,7,0,3,6,5,4), nrow=3)
	}else if (inout=='out'&& met=='tau'){
		flow<-matrix(rev(c(8,1,2,7,0,3,6,5,4)), nrow=3)
	}
		
	compare<-function(x){
		if(length(x[is.na(x)])==9){
			nmatch<-NA
		}else{
			comp<-x[!is.na(x)]==flow[!is.na(x)]
			nmatch<-length(comp[comp==TRUE])
		}
		nmatch[nmatch==0]<-1
		nmatch
	}
	
	r <- focal(flowdir, w=matrix(1,nrow=3,ncol=3), fun=compare, ...)
	r
	
}
