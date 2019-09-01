#' splash.grid
#'
#' Apply splash algorithm
#' @param   sw_in, lon
#' @param   tc, lon
#' @param   pn, lon
#' @param   elev, lon
#' @return a matrix xts type
#' @import Rcpp
#' @import raster 
#' @import parallel  
#' @import zoo
#' @keywords splash
#' @export
#' @examples
#' splash.grid()
splash.grid<-function(sw_in, tc, pn, elev, soil, outdir=getwd(),sim.control=list(par=TRUE, ncores=7,output.mode="monthly",inmem=FALSE), ...){
	#### IMPORT SOURCES ##########################################################
	# require(raster)
	# require(xts)
	# require(doSNOW)
	# require(zoo)
	rasterOptions(maxmemory=1e9, tmptime = 24, chunksize = 1e8,todisk = FALSE, overwrite=TRUE, tolerance = 0.5)
	
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
			cat("reaching steady state"); flush.console();
			beginCluster(sim.control$ncores, ...)
			eq<-spinup.grid(sw_in,tc,pn,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			cat(paste("solving","year",y[1])); flush.console();
			result.all<-run_one_year.grid(sw_in,tc,pn,eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			# endCluster()
			gc()
		}
		else if(length(y)>1){
			
			end<-cumsum(ny)
			start<-end+1
			result<-list()
			beginCluster(sim.control$ncores, ...)
			cat("reaching steady state"); flush.console();
			eq<-spinup.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]],elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			cat(paste("solving","year",y[1]))
			result[[1]]<-run_one_year.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]], eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],
				resolution,Au,sim.control$inmem,outdir=tmpdir)
			rm(eq)
			gc()
			###########################################################################
			# 3.2. loop through years
			###########################################################################
			# correct for leap years	inside c++ functions
			pb <- txtProgressBar(min=1,max = length(y), 3)
			for (i in 2:length(y)){
				cat(paste("solving","year",y[i])); flush.console();
				
				stidx<-i-1
				result[[i]]<-run_one_year.grid(sw_in[[start[stidx]:end[i]]], tc[[start[stidx]:end[i]]], pn[[start[stidx]:end[i]]],result[[stidx]]$wn,
					result[[stidx]]$snow,elev,lat,terraines,soil,y[i],resolution,Au,sim.control$inmem,outdir=tmpdir)
				setTxtProgressBar(pb,i)				
			}
			close(pb)
			# endCluster()
			
			gc()
			
		}
		
		
	}
	
	###########################################################################
	# 03. Start the calculations for Monthly inputs
	###########################################################################
	
	else if (abs(as.numeric(time.freq, units = "days"))>20){		
		
		if (length(y)==1){
			cat("reaching steady state"); flush.console();
			beginCluster(sim.control$ncores, ...)
			eq<-spinup.grid(sw_in,tc,pn,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			cat(paste("solving","year",y[1])); flush.console();
			result.all<-run_one_year.grid(sw_in,tc,pn,eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			# endCluster()
		}
		else if(length(y)>1){
			nm <- rep(12,length(y))
			end<-cumsum(nm)
			start<-end-11
			result<-list()
			beginCluster(sim.control$ncores, ...)
			cat("reaching steady state"); flush.console();
			eq<-spinup.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]],elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			cat(paste("solving","year",y[1])); flush.console();
			result[[1]]<-run_one_year.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]], eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],
				resolution,Au,sim.control$inmem,outdir=tmpdir)
			# endCluster()
			gc()
			# beginCluster(sim.control$ncores)
			###########################################################################
			# 3.2. loop through years
			###########################################################################
			# correct for leap years	inside c++ functions
			pb <- txtProgressBar(min=1,max = length(y),style = 3)
			
			for (i in 2:length(y)){
				cat(paste("solving","year",y[i])); flush.console();
				
				stidx<-i-1
				
				result[[i]]<-run_one_year.grid(sw_in[[start[i]:end[i]]], tc[[start[i]:end[i]]], pn[[start[i]:end[i]]],result[[stidx]]$wn,
					result[[stidx]]$snow,elev,lat,terraines,soil,y[i],resolution,Au,sim.control$inmem,outdir=tmpdir)
				
				setTxtProgressBar(pb,i)
				
			}
			close(pb)
			
			
			gc()
		}
		
	}		
	
	###########################################################################
	# 4. Building the raster stacks
	###########################################################################
	if(length(y)>1){
		flush.console();
		cat("building the grids"); flush.console();
		
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
		gc()
	}
	###########################################################################
	# 5. setting the time and saving to the disk
	###########################################################################
	settime<-function(x){
		x<-raster::setZ(x,ztime.days)
	}
	
	result.all<-lapply(result.all,settime)
	
	if (sim.control$output.mode=="monthly"){
		cat("aggregating to monthly"); flush.console();
		###########################################################################
		# 5. Define function to aggregate in parallel
		###########################################################################
		aggregate_par<-function(x,func='mean',ind.months=ztime.months,inmem=FALSE,varnam='wn',outdir=getwd()){
			###############################################################################################
			# 00. create array for results, fluxes: mm/day, storages (wn, snow): mm
			###############################################################################################
			y<-as.numeric(unique(format(ind.months,'%Y')))
			
			gc()
			nm <-length(ind.months)
			# actual soil moisture
			out<-brick(nrows=nrow(x), ncols=ncol(x), crs=crs(x), nl=nm)
			extent(out)<-extent(x)
			out<-setZ(out,ind.months)
			indmonth<-format(getZ(x),'%Y-%m')
			setwd(outdir)
			###############################################################################################
			# 01. set the clusters for parallel computing
			###############################################################################################	
			cl <- getCluster()
			on.exit( returnCluster() )
			nodes <- length(cl)
			bs <- blockSize(x, minblocks=nodes*5)
			parallel::clusterExport(cl, varlist=c('x','func','indmonth','bs'),envir=environment()) 
			pb <- pbCreate(bs$n)
			pb <- txtProgressBar(min=1,max = bs$n, style = 1)
			###############################################################################################
			# 02. create the functions to send to the workers, split the data in chunks
			###############################################################################################	
			clFun <- function(i) {
				
				x_block<-getValues(x,bs$row[i], bs$nrows[i])
				# do calculations
				if(func=='mean'){
					result<-t(sweep(rowsum(t(x_block),indmonth,na.rm = FALSE), 1, table(indmonth), "/") )
				}else{
					result<-t(rowsum(t(x_block),indmonth,na.rm = FALSE))
				}
				
				return(result)
			}
			###############################################################################################
			# 03. send tasks to the nodes
			###############################################################################################
			for (i in 1:nodes) {
				parallel:::sendCall(cl[[i]], clFun, i, tag=i)
			}
			###############################################################################################
			# 04. write to the disk on the go, or save to the ram
			###############################################################################################
			if(varnam=='wn'){
				longname='monthly soil moisture'
			}else if(varnam=='ro'){
				longname='monthly runoff'
			}else if(varnam=='pet'){
				longname='monthly potential evapotranspiration'
			}else if(varnam=='aet'){
				longname='monthly actual evapotranspiration'
			}else if(varnam=='snow'){
				longname='monthly snow water equivalent'
			}else if(varnam=='cond'){
				longname='monthly condensation'
			}else if(varnam=='bflow'){
				longname='monthly baseflow'
			}
			if(!inmem){
				out<-writeStart(out,filename=paste0(outdir,"/",y[1],"_",y[length(y)],".",varnam,".","nc"),format="CDF",overwrite=TRUE,varname=varnam, varunit="mm",
					longname=longname, xname="lon", yname="lat", zname="time", zunit=paste("months","since",paste0(y[1]-1,"-",12)))
				
				
			}else {
				matout <- matrix(ncol=nlayers(out), nrow=ncell(out))
				endind<-cumsum(bs$nrows*out@ncols)
				startind<-c(1,endind+1)    
			}
			###############################################################################################
			# 05. receive results from the nodes
			###############################################################################################	
			for (i in 1:bs$n) {
				
				d <- parallel:::recvOneData(cl)
				# error?
				if (! d$value$success) {
					stop('error!! check the data...')
				}
				# which block is this?
				b <- d$value$tag
				# cat('received block: ',b,'\n'); flush.console();
				if (!inmem) {
					out <- writeValues(out,d$value$value, bs$row[b])
					
				} else {
					
					matout[startind[b]:endind[b],] <- d$value$value
				}
				
				# need to send more data?
				ni <- nodes + i
				if (ni <= bs$n) {
					parallel:::sendCall(cl[[d$node]], clFun, ni, tag=ni)
				}
				setTxtProgressBar(pb,i)
			}
			###############################################################################################
			# 06. close connection with the files, or assign valueas to the raster objects
			###############################################################################################
			
			if (!inmem) {
				out <- writeStop(out)
				
			} else {
				
				out<-setValues(out,matout)
				
			}
			close(pb)
			gc()
			return(out)
			
		}
		###############################################################################################
		# 06. Do the aggregations one by one, using lapply causes a weird memory leak
		###############################################################################################
		
		result.all[[1]]<-aggregate_par(result.all[[1]],func='mean',ind.months=ztime.months,inmem=sim.control$inmem,varnam='wn',outdir=outdir)
		result.all[[2]]<-aggregate_par(result.all[[2]],func='sum',ind.months=ztime.months,inmem=sim.control$inmem,varnam='ro',outdir=outdir)
		result.all[[3]]<-aggregate_par(result.all[[3]],func='sum',ind.months=ztime.months,inmem=sim.control$inmem,varnam='pet',outdir=outdir)
		gc()
		result.all[[4]]<-aggregate_par(result.all[[4]],func='sum',ind.months=ztime.months,inmem=sim.control$inmem,varnam='aet',outdir=outdir)
		result.all[[5]]<-aggregate_par(result.all[[5]],func='mean',ind.months=ztime.months,inmem=sim.control$inmem,varnam='snow',outdir=outdir)
		result.all[[6]]<-aggregate_par(result.all[[6]],func='sum',ind.months=ztime.months,inmem=sim.control$inmem,varnam='cond',outdir=outdir)
		result.all[[7]]<-aggregate_par(result.all[[7]],func='sum',ind.months=ztime.months,inmem=sim.control$inmem,varnam='bflow',outdir=outdir)
		endCluster()
		
	}else{
		endCluster()
		gc()
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
	
	return(result.all)
	
}

spinup.grid<-function(sw_in, tc, pn, elev,lat, terraines,soil, y, resolution,  Au ,inmem=FALSE,outdir=getwd()){
	###############################################################################################
	# 00. create array for equilibrium soil moisture wneq and snow
	###############################################################################################
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	wneq<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(wneq)<-extent(elev)
	snoweq<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
	extent(snoweq)<-extent(elev)	
	setwd(outdir)
	if(inmem){
		sw_in<-readAll(sw_in)
		tc<-readAll(tc)
		pn<-readAll(pn)
		elev<-readAll(elev)
		lat<-readAll(lat)
		terraines<-readAll(terraines)
		soil<-readAll(soil)
		Au<-readAll(Au)
	}
	###############################################################################################
	# 01. set the clusters for parallel computing
	###############################################################################################	
	cl <- getCluster()
	on.exit( returnCluster() )
	nodes <- length(cl)
	message('Using cluster with ', nodes, ' nodes')
	bs <- blockSize(sw_in, minblocks=nodes*5)
	parallel::clusterExport(cl, c("sw_in","tc","pn","elev","lat","terraines",'soil','y','resolution','Au','bs'),envir=environment()) 
	pb <- pbCreate(bs$n)
	pb <- txtProgressBar(min=1,max = bs$n, style = 3)
	###############################################################################################
	# 02. create the functions to send to the workers, split the data in chunks
	###############################################################################################	
	clFun <- function(i) {
		swrow<- split(getValues(sw_in, bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		tcrow<-split(getValues(tc,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		pnrow<-split(getValues(pn,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		latrow<-getValues(lat,bs$row[i], bs$nrows[i])
		elevrow<-getValues(elev,bs$row[i], bs$nrows[i])
		sloprow<-getValues(terraines[[1]],bs$row[i], bs$nrows[i])
		asprow<-getValues(terraines[[2]],bs$row[i], bs$nrows[i])
		soilrow<-split(getValues(soil,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		Aurow<-getValues(Au,bs$row[i], bs$nrows[i])
		# do calculations
		wneqmat<-mapply(rspin_up,latrow,elevrow,swrow,tcrow,pnrow,sloprow,asprow,y,soilrow,Aurow,resolution)
		return(wneqmat)
	}
	###############################################################################################
	# 03. send tasks to the nodes
	###############################################################################################
	for (i in 1:nodes) {
		parallel:::sendCall(cl[[i]], clFun, i, tag=i)
	}
	###############################################################################################
	# 04. write to the disk on the go, or save to the ram
	###############################################################################################
	if(!inmem){
		wneq<-writeStart(wneq,filename="wneq.grd",overwrite=TRUE)
		snoweq<-writeStart(snoweq,filename="snoweq.grd",overwrite=TRUE)
	}else {
		matwneq <- matrix(ncol=nlayers(wneq), nrow=ncell(wneq))
		matswoeq<-matrix(ncol=nlayers(wneq), nrow=ncell(wneq))
		endind<-cumsum(bs$nrows*snoweq@ncols)
		startind<-c(1,endind+1)    
	}
	###############################################################################################
	# 05. receive results from the nodes
	###############################################################################################	
	for (i in 1:bs$n) {
		
		d <- parallel:::recvOneData(cl)
		# error?
		if (! d$value$success) {
			stop('cluster error')
		}
		# which block is this?
		b <- d$value$tag
		# cat('received block: ',b,'\n'); flush.console();
		if (!inmem) {
			wneq <- writeValues(wneq,do.call(rbind,d$value$value[1,]), bs$row[b])
			snoweq <- writeValues(snoweq, do.call(rbind,d$value$value[2,]), bs$row[b])
		} else {
			                
			matwneq[startind[b]:endind[b],] <- do.call(rbind,d$value$value[1,])
			matswoeq[startind[b]:endind[b],] <- do.call(rbind,d$value$value[2,])
		}
				
		# need to send more data?
		ni <- nodes + i
		if (ni <= bs$n) {
			parallel:::sendCall(cl[[d$node]], clFun, ni, tag=ni)
		}
		setTxtProgressBar(pb,i)
	}
	###############################################################################################
	# 06. close connection with the files, or assign valueas to the raster objects
	###############################################################################################
	
	if (!inmem) {
		wneq <- writeStop(wneq)
		snoweq <- writeStop(snoweq)
	} else {
		wneq<-setValues(wneq,matwneq)
		snoweq<-setValues(snoweq,matswoeq)
	}
	close(pb)
	gc()
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
	# make bricks, raster stacks are not working, result should be as stack, otherwise merging everithing wont work
	gc()
	setwd(outdir)
	if(inmem){
		sw_in<-readAll(sw_in)
		tc<-readAll(tc)
		pn<-readAll(pn)
		wn<-readAll(wn)
		snow<-readAll(snow)
		elev<-readAll(elev)
		lat<-readAll(lat)
		terraines<-readAll(terraines)
		soil<-readAll(soil)
		Au<-readAll(Au)
	}
	###############################################################################################
	# 01. set the clusters for parallel computing
	###############################################################################################	
	cl <- getCluster()
	on.exit( returnCluster() )
	nodes <- length(cl)
	bs <- blockSize(sw_in, minblocks=nodes*5)
	parallel::clusterExport(cl, varlist=c("sw_in","tc","pn",'wn','snow',"elev","lat","terraines",'soil','y','resolution','Au','bs'),envir=environment()) 
	pb <- pbCreate(bs$n)
	pb <- txtProgressBar(min=1,max = bs$n, style = 1)
	###############################################################################################
	# 02. create the functions to send to the workers, split the data in chunks
	###############################################################################################	
	clFun <- function(i) {
		swrow<- split(getValues(sw_in, bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		tcrow<-split(getValues(tc,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		pnrow<-split(getValues(pn,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		wnrow<-split(getValues(wn,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		snowrow<-split(getValues(snow,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		latrow<-getValues(lat,bs$row[i], bs$nrows[i])
		elevrow<-getValues(elev,bs$row[i], bs$nrows[i])
		sloprow<-getValues(terraines[[1]],bs$row[i], bs$nrows[i])
		asprow<-getValues(terraines[[2]],bs$row[i], bs$nrows[i])
		soilrow<-split(getValues(soil,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		Aurow<-getValues(Au,bs$row[i], bs$nrows[i])
		# do calculations
		yearlist<-mapply(run_one_year,latrow,elevrow,sloprow,asprow,swrow,tcrow,pnrow,wnrow,y,snowrow,soilrow,Aurow,resolution)
				
		return(yearlist)
	}
	###############################################################################################
	# 03. send tasks to the nodes
	###############################################################################################
	for (i in 1:nodes) {
		parallel:::sendCall(cl[[i]], clFun, i, tag=i)
	}
	###############################################################################################
	# 04. write to the disk on the go, or save to the ram
	###############################################################################################
	if(!inmem){
		sm<-writeStart(sm,filename=paste0("sm","_",y,".grd"),overwrite=TRUE)
		ro<-writeStart(ro,filename=paste0("ro","_",y,".grd"),overwrite=TRUE)
		pet<-writeStart(pet,filename=paste0("pet","_",y,".grd"),overwrite=TRUE)
		aet<-writeStart(aet,filename=paste0("aet","_",y,".grd"),overwrite=TRUE)
		swe<-writeStart(swe,filename=paste0("swe","_",y,".grd"),overwrite=TRUE)
		cond<-writeStart(cond,filename=paste0("cond","_",y,".grd"),overwrite=TRUE)
		bflow<-writeStart(bflow,filename=paste0("bflow","_",y,".grd"),overwrite=TRUE)
		
	}else {
		matsm <- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matro<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matpet <- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		mataet <- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matswe<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matcond <- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matbflow<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		endind<-cumsum(bs$nrows*sm@ncols)
		startind<-c(1,endind+1)    
	}
	###############################################################################################
	# 05. receive results from the nodes
	###############################################################################################	
	for (i in 1:bs$n) {
		
		d <- parallel:::recvOneData(cl)
		# error?
		if (! d$value$success) {
			stop('error!! check the data...')
		}
		# which block is this?
		b <- d$value$tag
		# cat('received block: ',b,'\n'); flush.console();
		if (!inmem) {
			sm <- writeValues(sm,do.call(rbind,d$value$value[1,]), bs$row[b])
			ro <- writeValues(ro, do.call(rbind,d$value$value[2,]), bs$row[b])
			pet <- writeValues(pet, do.call(rbind,d$value$value[3,]), bs$row[b])
			aet <- writeValues(aet, do.call(rbind,d$value$value[4,]), bs$row[b])
			swe <- writeValues(swe, do.call(rbind,d$value$value[5,]), bs$row[b])
			cond <- writeValues(cond, do.call(rbind,d$value$value[6,]), bs$row[b])
			bflow <- writeValues(bflow, do.call(rbind,d$value$value[7,]), bs$row[b])
			
		} else {
			
			matsm[startind[b]:endind[b],] <- do.call(rbind,d$value$value[1,])
			matro[startind[b]:endind[b],] <- do.call(rbind,d$value$value[2,])
			matpet[startind[b]:endind[b],] <- do.call(rbind,d$value$value[3,])
			mataet[startind[b]:endind[b],] <- do.call(rbind,d$value$value[4,])
			matswe[startind[b]:endind[b],] <- do.call(rbind,d$value$value[5,])
			matcond[startind[b]:endind[b],] <- do.call(rbind,d$value$value[6,])
			matbflow[startind[b]:endind[b],] <- do.call(rbind,d$value$value[7,])
		}
		
		# need to send more data?
		ni <- nodes + i
		if (ni <= bs$n) {
			parallel:::sendCall(cl[[d$node]], clFun, ni, tag=ni)
		}
		setTxtProgressBar(pb,i)
	}
	###############################################################################################
	# 06. close connection with the files, or assign valueas to the raster objects
	###############################################################################################
	
	if (!inmem) {
		sm <- writeStop(sm)
		ro <- writeStop(ro)
		pet <- writeStop(pet)
		aet <- writeStop(aet)
		swe <- writeStop(swe)
		cond <- writeStop(cond)
		bflow <- writeStop(bflow)	
	} else {
		# soil water content
		sm<-setValues(sm,matsm)
		# runoff
		ro<-setValues(ro,matro)
		# Potential evapotranspiration
		pet<-setValues(pet,matpet)
		# Actual evapotranspiration
		aet<-setValues(aet,mataet)
		# Snow water equivalent
		swe<-setValues(swe,matswe)
		# Condensation
		cond<-setValues(cond,matcond)
		# baseflow
		bflow<-setValues(bflow,matbflow)
		
	}
	close(pb)
	gc()
	return(list(wn=stack(sm),ro=stack(ro),pet=stack(pet),aet=stack(aet),snow=stack(swe),cond=stack(cond),bflow=stack(bflow)))
		
}
