#' upslope_area
#'
#' Computes the upslope area in m2 using topmodel for small areas or Taudem for large, corrections acording to latitude done by raster::area
#' @param   dem
#' @return a RasterStack with the contributing area per pixel, number of cells draining in, and number of cells draining out
#' @import topmodel
#' @import raster 
#' @keywords contributing area
#' @export
#' @examples
#' splash.grid()
upslope_area<-function(dem){

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
	upflow<-terrain(elev*-1,opt='flowdir')
	hg_sl<-tan(terrain(elev,opt='slope', unit='radians'))
	rm(elev)
	gc()
	# funct to extract the average upslope hydraulic gradient
	# calc_hg_path<-function(cell,fd,hg_local){
	# 	path<- flowPath(fd, cell)
	# 	if(is.null(path)){
	# 		return(NA)
	# 	}else{
	# 		hg<-.Internal(mean(extract(hg_local,path)))
	# 		return(hg)
	# 	}
	# }
	# 	
	# parallel:::clusterExport(getCluster(), c('calc_hg_path','upflow','hg_sl'),envir=environment()) 
	# 
	# hg_values<-parallel::clusterMap(cl = getCluster(), fun=calc_hg_path, cell=1:ncell(upflow),MoreArgs = list(fd=upflow,hg_local=hg_sl),SIMPLIFY = T)
	# 
	# hg<-raster(dem)
	# hg<-setValues(hg,hg_values)
	
	return(stack(areacatch,ncellin,ncellout))
	
}
upslope_areav2<-function(dem,type,tmpd){
	# returns upslope area in square meters
	# require(raster)
	# Set working directory to your location
	setwd(tmpd)
	resolution<-sqrt(cellStats(area(dem), stat='mean', na.rm=TRUE))*1000
	writeRaster(dem,"rawdem.tif",format="GTiff", overwrite=TRUE)
	
	if(type == 'MPIcluster'){
		# Pitremove
		system("pitremove -z rawdem.tif -fel dem_nopit.tif")
		
		# D8 flow directions
		system("d8flowdir -p dem_p.tif -sd8 dem_sd8.tif -fel dem_nopit.tif")
		# Contributing area
		system("aread8 -p dem_p.tif -ad8 dem_a_ac.tif -nc")
	}else{
		# Pitremove
		system("mpiexec pitremove -z rawdem.tif -fel dem_nopit.tif")
		
		# D8 flow directions
		system("mpiexec d8flowdir -p dem_p.tif -sd8 dem_sd8.tif -fel dem_nopit.tif")
		# Contributing area
		system("mpiexec aread8 -p dem_p.tif -ad8 dem_a_ac.tif -nc")
	}
			
	ups_ncell<-raster("dem_a_ac.tif")
	flowdir<-raster("dem_p.tif")
	hg_sl<-raster('dem_sd8.tif')
	hg_sl[hg_sl==0]<-0.001
	area_p_cell<-area(ups_ncell)
	ups_area<-overlay(ups_ncell, area_p_cell, fun=function(x,y){(x*y*1e6)},filename="areacatch.grd",overwrite=TRUE)
	ncellin<-ncellflow(flowdir,inout='in',met='tau',filename="ncellin.grd",overwrite=TRUE)
	ncellout<-ncellflow(flowdir,inout='out',met='tau',filename="ncellout.grd",overwrite=TRUE)
	# upflow<-terrain(raster("dem_nopit.tif")*-1,opt='flowdir')
	# # funct to extract the average upslope hydraulic gradient
	# calc_hg_path<-function(cell,fd,hg_local){
	# 	path<- flowPath(fd, cell)
	# 	if(is.null(path)){
	# 		return(NA)
	# 	}else{
	# 		hg<-.Internal(mean(extract(hg_local,path)))
	# 		return(hg)
	# 	}
	# }
	# 
	# parallel:::clusterExport(getCluster(), c('calc_hg_path','upflow','hg_sl'),envir=environment()) 
	# 
	# hg_values<-parallel::clusterMap(cl = getCluster(), fun=calc_hg_path, cell=1:ncell(upflow),MoreArgs = list(fd=upflow,hg_local=hg_sl),SIMPLIFY = T)
	# 
	# hg<-raster(dem)
	# hg<-setValues(hg,hg_values)
	# hg<-log((ups_area)/hg_sl)
		
	return(stack(ups_area,ncellin,ncellout))
	
	
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
