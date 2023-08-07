#' upslope_area
#'
#' Computes n cells draining in/out and the contributing area per pixel in m2 using topmodel for small areas or Taudem for large, corrections acording to latitude done by raster::area
#' @param   dem
#' @return a RasterStack with the contributing area per pixel, number of cells draining in, and number of cells draining out
#' @import topmodel
#' @import raster 
#' @keywords contributing area
#' @export
#' @examples
#' splash.grid()
upslope_area<-function(dem,resolution){

	resolution<-cellStats(resolution, stat='mean', na.rm=TRUE)
	#1.0.fill the sinks
	elev<-sinkfill(raster::as.matrix(dem),resolution,1)
	elev<-raster(elev)
	extent(elev)<-extent(dem)
	crs(elev)<-crs(dem)
	#2. calculate contributing area
	areacatch<-topidx(raster::as.matrix(elev), resolution= resolution)$area
	areacatch<-raster(areacatch)
	extent(areacatch)<-extent(elev)
	crs(areacatch)<-crs(elev)
	areacatch<-writeRaster(areacatch,"areacatch.grd",overwrite=TRUE)
	#3. calculate flow direction and n cells draining in and out
	flowdir<-terrain(elev,opt='flowdir')
	ncellin<-ncellflow(flowdir,inout='in',met='top',filename="ncellin.grd",overwrite=TRUE)
	ncellout<-ncellflow(flowdir,inout='out',met='top',filename="ncellout.grd",overwrite=TRUE)
	
	rm(elev)
	gc()
	###################################################################################
	#experimental, upslope average hydraulic gradient, flow paths and flow distance, too slow for now
	###################################################################################
	#upflow<-terrain(elev*-1,opt='flowdir')
	#hg_sl<-tan(terrain(elev,opt='slope', unit='radians'))
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
	
	return(raster::stack(areacatch,ncellin,ncellout))
	
}
upslope_areav2<-function(dem,type,tmpd){
	# returns upslope area in square meters
	# require(raster)
	# Set working directory to your location
	setwd(tmpd)
	writeRaster(dem,"rawdem.tif",datatype='INT2S',format="GTiff", overwrite=TRUE)
		
	if(type != 'SOCKcluster'){
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
	#2. compute contributing area
	ups_ncell<-raster("dem_a_ac.tif")
	area_p_cell<-area(ups_ncell)
	
	###fix coastal regions with flat slopes
	fix_coast<-function(elev,ups_ncell){
		ups_ncell<-ifelse(!is.na(elev) & is.na(ups_ncell), 1.0, ups_ncell)
		ups_ncell
	}
	
	ups_ncell<-overlay(dem,ups_ncell,fun=fix_coast,filename="ups_ncell.grd",overwrite=TRUE)
	### calculate upslope area
	ups_area<-overlay(ups_ncell, area_p_cell, fun=function(x,y){(x*y*1e6)},filename="areacatch.grd",overwrite=TRUE)
	#3. calculate flow direction and n cells draining in and out
	flowdir<-raster("dem_p.tif")
		
	if(ncell(dem)>1e7){
		doParallel::registerDoParallel(getCluster())
		ncellin<-ncellflow_par(flowdir,inout='in',met='tau',filename="ncellin.grd",overwrite=TRUE)
		ncellout<-ncellflow_par(flowdir,inout='out',met='tau',filename="ncellout.grd",overwrite=TRUE)
		gc()
	}else{
		ncellin<-ncellflow(flowdir,inout='in',met='tau',filename="ncellin.grd",overwrite=TRUE)
		ncellout<-ncellflow(flowdir,inout='out',met='tau',filename="ncellout.grd",overwrite=TRUE)
	}
	
	###################################################################################
	#experimental, upslope average hydraulic gradient, flow paths and flow distance, too slow for now
	###################################################################################
	# hg_sl<-raster('dem_sd8.tif')
	# hg_sl[hg_sl==0]<-0.001
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
		
	return(raster::stack(ups_area,ncellin,ncellout))
	
	
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


ncellflow_par<-function(flowdir,inout='in',met='top', ...){
	
	
	compare<-function(x,inout,met){
		if(inout=='in'&& met=='top'){
			flow<-matrix(c(2,1,128,4,0,64,8,16,32), nrow=3)
		}else if (inout=='out'&& met=='top'){
			flow<-matrix(rev(c(2,1,128,4,0,64,8,16,32)), nrow=3)
		}else if (inout=='in'&& met=='tau'){
			flow<-matrix(c(8,1,2,7,0,3,6,5,4), nrow=3)
		}else if (inout=='out'&& met=='tau'){
			flow<-matrix(rev(c(8,1,2,7,0,3,6,5,4)), nrow=3)
		}
		if(length(x[is.na(x)])==9){
			nmatch<-NA
		}else{
			comp<-x[!is.na(x)]==flow[!is.na(x)]
			nmatch<-length(comp[comp==TRUE])
		}
		nmatch[nmatch==0]<-1
		nmatch
	}
	
	# r <- focal(flowdir, w=matrix(1,nrow=3,ncol=3), fun=compare, ...)
	# sfQuickInit(cpus=6)
	# registerDoParallel(get)
	r <- spatial.tools::focal_hpc(x=flowdir,fun=compare,args=list(inout=inout,met=met),window_dims=c(3,3))
	r
	
}












