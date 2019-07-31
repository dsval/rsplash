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


run_one_year.grid<-function(sw_in, tc, pn,wn,snow ,elev,lat, terraines,soil, y, resolution,  Au,inmem=FALSE,outdir=getwd()){
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
		setwd(outdir)
		sm<-writeStart(sm,filename=paste0("sm","_",y,".grd"),overwrite=TRUE)
		ro<-writeStart(ro,filename=paste0("ro","_",y,".grd"),overwrite=TRUE)
		pet<-writeStart(pet,filename=paste0("pet","_",y,".grd"),overwrite=TRUE)
		aet<-writeStart(aet,filename=paste0("aet","_",y,".grd"),overwrite=TRUE)
		swe<-writeStart(swe,filename=paste0("swe","_",y,".grd"),overwrite=TRUE)
		cond<-writeStart(cond,filename=paste0("cond","_",y,".grd"),overwrite=TRUE)
		bflow<-writeStart(bflow,filename=paste0("bflow","_",y,".grd"),overwrite=TRUE)
		
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
	# soil moisture and swe to the memory, otherwise painfully slow
	sm<-readAll(sm)
	swe<-readAll(swe)
	
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
	tmpdir<-dirname(rasterTmpFile())
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
			eq<-spinup.grid(sw_in,tc,pn,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			cat(paste("solving","year",y[1]))
			result.all<-run_one_year.grid(sw_in,tc,pn,eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			stopCluster(cl)
		}
		else if(length(y)>1){
			
			end<-cumsum(ny)
			start<-end+1
			result<-list()
			cl <- parallel::makeCluster(sim.control$ncores, ...)
			doSNOW::registerDoSNOW(cl)
			cat("reaching steady state")
			eq<-spinup.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]],elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			cat(paste("solving","year",y[1]))
			result[[1]]<-run_one_year.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]], eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],
										resolution,Au,sim.control$inmem,outdir=tmpdir)
			
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
					result[[stidx]]$snow,elev,lat,terraines,soil,y[i],resolution,Au,sim.control$inmem,outdir=tmpdir)				
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
			eq<-spinup.grid(sw_in,tc,pn,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			cat(paste("solving","year",y[1]))
			result.all<-run_one_year.grid(sw_in,tc,pn,eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
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
			eq<-spinup.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]],elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			cat(paste("solving","year",y[1]))
			result[[1]]<-run_one_year.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]], eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],
										resolution,Au,sim.control$inmem,outdir=tmpdir)
			
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
					result[[stidx]]$snow,elev,lat,terraines,soil,y[i],resolution,Au,sim.control$inmem,outdir=tmpdir)
				
				
			}
			close(pb)
			stopCluster(cl)
			
			gc()
		}
		
	}		
	
	###########################################################################
	# 4. Building the raster stacks
	###########################################################################
	if(length(y)>1){
		
		cat("building the grids")
		
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
	
	
	
	if (sim.control$output.mode=="monthly"){
		cat("aggregating to monthly")
		###########################################################################
		# 5. Aggregate monthly using parallel processing
		###########################################################################
		# define monthly index for tapply
		indmonth<-format(ztime.days,'%Y-%m')
		# define functions to aggregate, not working with ratertempfile() use timestamp instead
		month_mean<-function(x){
			
			x<-calc(x=x, fun=function(x)tapply(x,indmonth,FUN= mean),filename=paste0(format(Sys.time(), "%Y%m%d_%H%M%S_"), "monthagg.grd"))
		}
		month_sum<-function(x){
			
			x<-calc(x=x, fun=function(x)tapply(x,indmonth,FUN= sum),filename=paste0(format(Sys.time(), "%Y%m%d_%H%M%S_"), "monthagg.grd"))
		}
		# serial slow for bigrasters
		# start.time<-Sys.time()
		# result.all[c(1,5)]<-lapply(result.all[c(1,5)], month_mean)
		# result.all[c(2,3,4,6,7)]<-lapply(result.all[c(2,3,4,6,7)], month_sum)
		# end.time<-Sys.time()
		# end.time-start.time
		# start.time<-Sys.time()
		cl <- parallel::makeCluster(sim.control$ncores, ...)
		doSNOW::registerDoSNOW(cl)
		snow::clusterEvalQ(cl, library("raster"))
		snow::clusterExport(cl, list=c("result.all","month_sum","month_mean","indmonth"),envir=environment()) 
		result.all[c(1,5)]<-snow::parLapply(cl,result.all[c(1,5)],month_mean)
		# memory leak, close cluster and cleanmem
		stopCluster(cl)
		gc()
		cl <- parallel::makeCluster(sim.control$ncores, ...)
		doSNOW::registerDoSNOW(cl)
		snow::clusterEvalQ(cl, library("raster"))
		snow::clusterExport(cl, list=c("result.all","month_sum","month_mean","indmonth"),envir=environment()) 
		result.all[c(2,3,4,6,7)]<-snow::parLapply(cl,result.all[c(2,3,4,6,7)],month_sum)
		stopCluster(cl)
		gc()
		# setting ztime
		settime<-function(x){
			x<-raster::setZ(x,ztime.months)
		}
		
		result.all<-lapply(result.all,settime)
		# end.time<-Sys.time()
		# end.time-start.time
		
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
	on.exit(stopCluster(cl))	
	return(result.all)
	
}
