#' Simple process-led algorithms for simulating habitats (SPLASH v.2.0)
#' 
#' R/C++ implementation of the SPLASH v.2.0 algorithm (Davis et al., 2017; Sandoval et al., in prep.).
#' 
#' @param   sw_in Incoming shortwave solar radiation (W m-2), Raster* object of monthly or daily averages with z time dimension.
#' @param   tc Air temperature (°C), same dimensions as sw_in
#' @param   pn Precipitation (mm), same dimensions as sw_in
#' @param   elev Elevation (m.a.s.l)
#' @param   soil Raster* object with the layers organized as sand(perc),clay(perc),organic matter(perc),coarse-fragments-fraction(perc), bulk density(g cm-3) and depth(m)
#' @param   outdir (optional) directory path where the results will be saved, working directory by default
#' @param   tmpdir (optional) directory path where the temporary files will be saved, default temporary directory by default
#' @param   sim.control (optional) list including options to control the output: output.mode="monthly" by default, "daily" also available, inmem=FALSE by default write all the results to the disk by chunks to save RAM, sacrificing speed.
#' @return a list of rasterBricks objects with z time dimension, all of them saved to outdir as netcdf files:
#' \itemize{
#'         \item \eqn{W_n}: Soil water content (mm) within the first 2 m of depth.
#'         \item \eqn{ro}: Runoff (mm d-1).
#'         \item \eqn{pet}: Potential evapotranspiration (mm d-1).
#'         \item \eqn{aet}: Actual evapotranspiration (mm d-1).
#'         \item \eqn{snow}: Snow water equivalent (mm).
#'         \item \eqn{bflow}: Lateral flow (mm d-1).
#' }
#' @import Rcpp 
#' @import xts
#' @keywords splash, evapotranspiration, soil moisture
#' @export
#' @examples
#' splash.grid(sw_in=200, tc=15, pn=10, lat=44,elev=1800,slop=10,asp=270,soil_data=c(sand=44,clay=2,OM=6,fgravel=12))
splash.grid<-function(sw_in, tc, pn, elev, soil, outdir=getwd(),tmpdir=dirname(rasterTmpFile()),sim.control=list(output.mode="monthly",inmem=FALSE)){
	###########################################################################
	# 00. Check if parallel computation is required by the user and if the dimensions of the raster objects match
	###########################################################################
	on.exit(endCluster())
	clcheck<-try(getCluster(), silent=TRUE)
	if(class(clcheck)=="try-error"){
		# If no cluster is initialized, assume only one core will do the calculations, beginCluster(1) saved me the time of coding serial versions of the functions
		beginCluster(1,'SOCK')
		message('Only using one core, use first beginCluster() if you want to run splash in parallel!!')
		
	}
	rasterOptions(tolerance = 0.5)
	compareRaster(sw_in, tc, pn, elev, soil,extent=TRUE, crs=TRUE, res=TRUE, orig=FALSE,rotation=FALSE, values=FALSE, stopiffalse=FALSE, showwarning=TRUE)	
	type=class(clcheck)[1]
	###########################################################################
	# 01. Calculate spatial distributed variables
	###########################################################################
		
	setwd(tmpdir)
		
	#### function to get the latitudes from big rasters i.e 1km res global extent
	getlatitude <- function(x, filename, ...) {
		##create array for the results
		out <-raster(x)
		# get the index of the blocks, every block has n rows, bigger the minblocks, smaller the chunk of rows
		bs <- blockSize(x, minblocks=200)
		pb <- pbCreate(bs$n)
		pb <- txtProgressBar(min=1,max = bs$n, style = 1)
		# start writing the outputfile
		out <- writeStart(out, filename, overwrite=TRUE)
		# rsqv<-function(x,y){summary(lm(y~x,na.action=))$r.squared}
		for (i in 1:bs$n) {
			# i=178
			# xmat <- getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
			# ymat <- getValues(y, row=bs$row[i], nrows=bs$nrows[i] )
			xncells<-cellFromRow(x,bs$row[i]:(bs$row[i]+ bs$nrows[i]-1))
			xmat<-getValues(x,bs$row[i], bs$nrows[i])
			xydata<-xyFromCell(x, xncells)
			
			# write the chunk of results, bs$row[i] is putting the results in the correct rows
			out <- writeValues(out, xydata[,2], bs$row[i])
			setTxtProgressBar(pb,i)
		}
		out <- writeStop(out)
		close(pb)
		return(out)
	}
	# 1.1 calculate upslope area in m2, *all the rasters will be saved to the disk by default
	if (ncell(elev)>1e7|res(elev)[1]>0.1){
		cat("computing contributing areas and terrain features, it might take a while...","\n")
		start.time<-Sys.time()
		Au<-upslope_areav2(elev, type,tmpdir)
		# 1.2.2 get latitude from the high res dem
		lat<-getlatitude(elev, filename='lat.grd')
		system("gdaldem slope -s 111120 -co BIGTIFF=YES rawdem.tif slope_deg.tif")
		system("gdaldem aspect -zero_for_flat -co BIGTIFF=YES rawdem.tif aspect_deg.tif")
		terraines<-raster::stack(list(slope='slope_deg.tif',aspect='aspect_deg.tif'))
		# get resolution in m2
		resolution<-area(elev,filename='area.grd',overwrite=T)
		resolution<-calc(resolution,function(x){sqrt(x)*1000})
		gc()
		cat("... it took","\n")
		end.time<-Sys.time()
		end.time-start.time
		cat("...","\n")
		
	}else{
		# get resolution in m2
		cat("computing contributing areas and terrain features","\n")
		resolution<-sqrt(area(elev))*1000
		Au<-upslope_area(elev,resolution)
		# 1.2 get latitudes
		lat<-elev*0
		lat.data<-rasterToPoints(elev)
		lat[!is.na(lat)]<-lat.data[,2]
		rm(lat.data)
		terraines<-terrain(elev, opt=c('slope', 'aspect'), unit='degrees')
		cat("done!","\n")
		gc()
	}
		
	
	#lat<-writeRaster(lat,filename="lat.grd", overwrite=TRUE)
	
	# 1.3  calculate slope and aspect
	
	
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
			cat("reaching steady state","\n")
			# beginCluster(sim.control$ncores,type=type)
			eq<-spinup.grid(sw_in,tc,pn,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			cat(paste("solving","year",y[1]),"\n")
			result.all<-run_one_year.grid(sw_in,tc,pn,eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],resolution,Au,eq$bfloweq,eq$tdraineq,sim.control$inmem,outdir=tmpdir)
			result.all$tdrain<-NULL
			# endCluster()
			gc()
		}
		else if(length(y)>1){
			
			end<-cumsum(ny)
			start<-end+1
			result<-list()
			# beginCluster(sim.control$ncores, type=type)
			cat("reaching steady state","\n")
			eq<-spinup.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]],elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			cat(paste("solving","year",y[1]),"\n")
			result[[1]]<-run_one_year.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]], eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],
				resolution,Au,eq$bfloweq,eq$tdraineq,sim.control$inmem,outdir=tmpdir)
			rm(eq)
			gc()
			###########################################################################
			# 3.2. loop through years
			###########################################################################
			# correct for leap years	inside c++ functions
			pb <- txtProgressBar(min=1,max = length(y), 3)
			for (i in 2:length(y)){
				cat(paste("solving","year",y[i]),"\n")
				
				stidx<-i-1
				result[[i]]<-run_one_year.grid(sw_in[[start[stidx]:end[i]]], tc[[start[stidx]:end[i]]], pn[[start[stidx]:end[i]]],result[[stidx]]$wn,
					result[[stidx]]$snow,elev,lat,terraines,soil,y[i],resolution,Au,result[[stidx]]$q_in,result[[stidx]]$tdrain,sim.control$inmem,outdir=tmpdir)
				setTxtProgressBar(pb,i)
				cat("...","\n")				
			}
			close(pb)
			# endCluster()
			
			gc()
			
		}
		
		
	}
	
	###########################################################################
	# 03. Start the calculations if the inputs are monthly
	###########################################################################
	
	else if (abs(as.numeric(time.freq, units = "days"))>20){		
		
		if (length(y)==1){
			cat("reaching steady state","\n")
			# beginCluster(sim.control$ncores, type=type)
			eq<-spinup.grid(sw_in,tc,pn,elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			cat(paste("solving","year",y[1]),"\n")
			result.all<-run_one_year.grid(sw_in,tc,pn,eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],resolution,Au,eq$bfloweq,eq$tdraineq, sim.control$inmem,outdir=tmpdir)
			result.all$tdrain<-NULL
			# endCluster()
		}
		else if(length(y)>1){
			nm <- rep(12,length(y))
			end<-cumsum(nm)
			start<-end-11
			result<-list()
			# beginCluster(sim.control$ncores, type=type)
			cat("reaching steady state","\n")
			eq<-spinup.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]],elev,lat,terraines,soil,y[1],resolution,Au,sim.control$inmem,outdir=tmpdir)
			cat(paste("solving","year",y[1]),"\n")
			result[[1]]<-run_one_year.grid(sw_in[[1:end[1]]], tc[[1:end[1]]], pn[[1:end[1]]], eq$wneq,eq$snoweq,elev,lat,terraines,soil,y[1],
				resolution,Au,eq$bfloweq,eq$tdraineq,sim.control$inmem,outdir=tmpdir)
			# endCluster()
			gc()
			# beginCluster(sim.control$ncores)
			###########################################################################
			# 3.2. loop through years
			###########################################################################
			# correct for leap years	inside c++ functions
			pb <- txtProgressBar(min=1,max = length(y),style = 3)
			
			for (i in 2:length(y)){
				cat(paste("solving","year",y[i]),"\n")
				
				stidx<-i-1
				
				result[[i]]<-run_one_year.grid(sw_in[[start[i]:end[i]]], tc[[start[i]:end[i]]], pn[[start[i]:end[i]]],result[[stidx]]$wn,
					result[[stidx]]$snow,elev,lat,terraines,soil,y[i],resolution,Au,result[[stidx]]$q_in,result[[stidx]]$tdrain,sim.control$inmem,outdir=tmpdir)
				
				setTxtProgressBar(pb,i)
				cat("...","\n")
				
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
		cat("building the grids","\n")
		
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
		# result.all$cond<-result[[1]]$cond
		# result.all$cond@layers<-purrr::flatten(lapply(result, function(x) {as.list(x$cond)}))
		result.all$bflow<-result[[1]]$bflow
		result.all$bflow@layers<-purrr::flatten(lapply(result, function(x) {as.list(x$bflow)}))
		# result.all$tdrain<-result[[1]]$tdrain
		# result.all$tdrain@layers<-purrr::flatten(lapply(result, function(x) {as.list(x$tdrain)}))
		
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
		cat("aggregating to monthly","\n")
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
			# save the empty array
			
			out<-brick(x[[1]], values=FALSE, nl=nm, filename=paste0(outdir,"/",y[1],"_",y[length(y)],".",varnam,".","nc"),format="CDF",overwrite=TRUE,
				varnam= varnam, xname="lon", yname="lat", zname="time")
				
			out<-setZ(out,ind.months)
			indmonth<-format(getZ(x),'%Y-%m')
			setwd(outdir)
			###############################################################################################
			# 01. set the clusters for parallel computing
			###############################################################################################	
			cl <- getCluster()
			on.exit( returnCluster() )
			nodes <- length(cl)
			bs <- blockSize(x, minblocks=nodes)
			parallel:::clusterExport(cl, c('x','func','indmonth','bs'),envir=environment()) 
			pb <- pbCreate(bs$n)
			pb <- txtProgressBar(min=1,max = max(bs$n,2), style = 1)
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
			}else if(varnam=='tdrain'){
				longname='days to drain'
			}
			if(!inmem){
				out<-writeStart(out,filename=paste0(outdir,"/",y[1],"_",y[length(y)],".",varnam,".","nc"),format="CDF",overwrite=TRUE,varname=varnam, varunit="mm",	longname=longname, xname="lon", yname="lat", zname="time")
				
				
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
		
		result.all$wn<-aggregate_par(result.all$wn,func='mean',ind.months=ztime.months,inmem=sim.control$inmem,varnam='wn',outdir=outdir)
		result.all$ro<-aggregate_par(result.all$ro,func='sum',ind.months=ztime.months,inmem=sim.control$inmem,varnam='ro',outdir=outdir)
		result.all$pet<-aggregate_par(result.all$pet,func='sum',ind.months=ztime.months,inmem=sim.control$inmem,varnam='pet',outdir=outdir)
		gc()
		result.all$aet<-aggregate_par(result.all$aet,func='sum',ind.months=ztime.months,inmem=sim.control$inmem,varnam='aet',outdir=outdir)
		result.all$snow<-aggregate_par(result.all$snow,func='mean',ind.months=ztime.months,inmem=sim.control$inmem,varnam='snow',outdir=outdir)
		# result.all[[6]]<-aggregate_par(result.all[[6]],func='sum',ind.months=ztime.months,inmem=sim.control$inmem,varnam='cond',outdir=outdir)
		result.all$bflow<-aggregate_par(result.all$bflow,func='sum',ind.months=ztime.months,inmem=sim.control$inmem,varnam='bflow',outdir=outdir)
		# result.all[[8]]<-aggregate_par(result.all[[8]],func='mean',ind.months=ztime.months,inmem=sim.control$inmem,varnam='tdrain',outdir=outdir)
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
		# result.all$cond<-writeRaster(result.all$cond,paste0(outdir,"/",y[1],"_",y[length(y)],".","cond",".","nc"),format="CDF",overwrite=TRUE,varname="cond", 
		# 	varunit="mm/day", longname="daily condensation water", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31)))
		result.all$bflow<-writeRaster(result.all$bflow,paste0(outdir,"/",y[1],"_",y[length(y)],".","bflow",".","nc"),format="CDF",overwrite=TRUE,varname="bflow", 
			varunit="mm/day", longname="daily baseflow", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31)))
		# result.all$tdrain<-writeRaster(result.all$tdrain,paste0(outdir,"/",y[1],"_",y[length(y)],".","tdrain",".","nc"),format="CDF",overwrite=TRUE,varname="tdrain", 
		# 	varunit="day", longname="days to drain", xname="lon", yname="lat", zname="time", zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31)))
		
	}
	
	gc()	
	
	return(result.all)
	
}

spinup.grid<-function(sw_in, tc, pn, elev,lat, terraines,soil, y, resolution,  Au ,inmem=FALSE,outdir=getwd()){
	###############################################################################################
	# 00. create array for equilibrium soil moisture wneq and snow
	###############################################################################################
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	setwd(outdir)
	if(!inmem){
		# actual soil moisture
		wneq<-brick(elev, values=FALSE, nl=ny,filename="wneq.grd",overwrite=TRUE,overwrite=TRUE) 
		snoweq<-brick(elev, values=FALSE, nl=ny,filename="snoweq.grd",overwrite=TRUE,overwrite=TRUE) 
		bfloweq<-brick(elev, values=FALSE, nl=ny,filename="bfloweq.grd",overwrite=TRUE,overwrite=TRUE) 
		tdraineq<-brick(elev, values=FALSE, nl=ny,filename="tdraineq.grd",overwrite=TRUE,overwrite=TRUE) 
		ro<-brick(elev, values=FALSE, nl=ny,filename="roeq.grd",overwrite=TRUE,overwrite=TRUE) 
	}else{
		
		wneq<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
		extent(wneq)<-extent(elev)
		snoweq<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
		extent(snoweq)<-extent(elev)
		bfloweq<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
		extent(bfloweq)<-extent(elev)
		tdraineq<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
		extent(tdraineq)<-extent(elev)
		# runoff
		ro<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
		extent(ro)<-extent(elev)
		#send the inputs to the RAM
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
	bs <- blockSize(sw_in, minblocks=nodes)
	parallel:::clusterExport(cl, c("sw_in","tc","pn","elev","lat","terraines",'soil','y','resolution','Au','bs'),envir=environment()) 
	pb <- pbCreate(bs$n)
	pb <- txtProgressBar(min=1,max = max(bs$n,2), style = 3)
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
		Aurow<-split(getValues(Au,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		resrow<-getValues(resolution,bs$row[i], bs$nrows[i])
		# do calculations
		wneqmat<-mapply(rspin_up,latrow,elevrow,swrow,tcrow,pnrow,sloprow,asprow,y,soilrow,Aurow,resrow)
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
		bfloweq<-writeStart(bfloweq,filename="bfloweq.grd",overwrite=TRUE)
		tdraineq<-writeStart(tdraineq,filename="tdraineq.grd",overwrite=TRUE)
		ro<-writeStart(ro,filename="roeq.grd",overwrite=TRUE)
	}else {
		matwneq <- matrix(ncol=nlayers(wneq), nrow=ncell(wneq))
		matswoeq<-matrix(ncol=nlayers(wneq), nrow=ncell(wneq))
		matbfeq<-matrix(ncol=nlayers(wneq), nrow=ncell(wneq))
		mattdeq<-matrix(ncol=nlayers(wneq), nrow=ncell(wneq))
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
			bfloweq <- writeValues(bfloweq, do.call(rbind,d$value$value[3,]), bs$row[b])
			tdraineq <- writeValues(tdraineq, do.call(rbind,d$value$value[4,]), bs$row[b])
			ro <- writeValues(ro, do.call(rbind,d$value$value[5,]), bs$row[b])
			
		} else {
			                
			matwneq[startind[b]:endind[b],] <- do.call(rbind,d$value$value[1,])
			matswoeq[startind[b]:endind[b],] <- do.call(rbind,d$value$value[2,])
			matbfeq[startind[b]:endind[b],] <- do.call(rbind,d$value$value[3,])
			mattdeq[startind[b]:endind[b],] <- do.call(rbind,d$value$value[4,])
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
		bfloweq <- writeStop(bfloweq)
		tdraineq <- writeStop(tdraineq)
		ro <- writeStop(ro)
	} else {
		wneq<-setValues(wneq,matwneq)
		snoweq<-setValues(snoweq,matswoeq)
		bfloweq<-setValues(bfloweq,matbfeq)
		tdraineq<-setValues(tdraineq,mattdeq)
	}
	close(pb)
	gc()
	return(list(wneq=wneq,snoweq=snoweq,bfloweq=bfloweq,tdraineq=tdraineq,ro=ro))
}


run_one_year.grid<-function(sw_in, tc, pn,wn,snow ,elev,lat, terraines,soil, y, resolution,  Au,bf_in, tdin,inmem=FALSE,outdir=getwd()){
	###############################################################################################
	# 00. create array for results, fluxes: mm/day, storages (wn, snow): mm 
	# *make bricks, raster stacks are not working, result should stacks, otherwise merging everything wont work
	###############################################################################################
	setwd(outdir)
	ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
	if(!inmem){
		# actual soil moisture
		sm<-brick(elev, values=FALSE, nl=ny,filename=paste0("sm","_",y,".grd"),overwrite=TRUE) 
		# runoff
		ro<-brick(elev, values=FALSE, nl=ny,filename=paste0("ro","_",y,".grd"),overwrite=TRUE) 
		# Potential evapotranspiration
		pet<-brick(elev, values=FALSE, nl=ny,filename=paste0("pet","_",y,".grd"),overwrite=TRUE) 
		# Actual evapotranspiration
		aet<-brick(elev, values=FALSE, nl=ny,filename=paste0("aet","_",y,".grd"),overwrite=TRUE) 
		# Snow water equivalent
		swe<-brick(elev, values=FALSE, nl=ny,filename=paste0("swe","_",y,".grd"),overwrite=TRUE) 
		# Condensation
		# cond<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
		# extent(cond)<-extent(elev)
		# baseflow
		bflow<-brick(elev, values=FALSE, nl=ny,filename=paste0("bflow","_",y,".grd"),overwrite=TRUE)
		# tdrain
		tdrain<-brick(elev, values=FALSE, nl=ny,filename=paste0("tdrain","_",y,".grd"),overwrite=TRUE)
		# q_in
		q_in<-brick(elev, values=FALSE, nl=ny,filename=paste0("q_in","_",y,".grd"),overwrite=TRUE)
	}else{
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
		# cond<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
		# extent(cond)<-extent(elev)
		# baseflow
		bflow<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=ny)
		extent(bflow)<-extent(elev)
		tdrain<-bflow
		q_in<-bflow
		#move inputs to the memory
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
		bf_in<-readAll(bf_in)
		tdin<-readAll(tdin)
	
	}
			
	###############################################################################################
	# 01. set the clusters for parallel computing
	###############################################################################################	
	cl <- getCluster()
	on.exit( returnCluster() )
	nodes <- length(cl)
	bs <- blockSize(sw_in, minblocks=nodes)
	parallel:::clusterExport(cl, c("sw_in","tc","pn",'wn','snow',"elev","lat","terraines",'soil','y','resolution','Au','bf_in', 'tdin','bs'),envir=environment()) 
	pb <- pbCreate(bs$n)
	pb <- txtProgressBar(min=1,max = max(bs$n,2), style = 1)
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
		Aurow<-split(getValues(Au,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		qinrow<-split(getValues(bf_in,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		tdrow<-split(getValues(tdin,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		resrow<-getValues(resolution,bs$row[i], bs$nrows[i])
		# do calculations
		yearlist<-mapply(run_one_year,latrow,elevrow,sloprow,asprow,swrow,tcrow,pnrow,wnrow,y,snowrow,soilrow,Aurow,resrow,qinrow,tdrow)
				
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
		# cond<-writeStart(cond,filename=paste0("cond","_",y,".grd"),overwrite=TRUE)
		bflow<-writeStart(bflow,filename=paste0("bflow","_",y,".grd"),overwrite=TRUE)
		tdrain<-writeStart(tdrain,filename=paste0("tdrain","_",y,".grd"),overwrite=TRUE)
		q_in<-writeStart(q_in,filename=paste0("q_in","_",y,".grd"),overwrite=TRUE)
		
	}else {
		matsm <- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matro<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matpet <- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		mataet <- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matswe<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		# matcond <- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matbflow<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		mattdrain<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matq_in<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
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
			# cond <- writeValues(cond, do.call(rbind,d$value$value[6,]), bs$row[b])
			bflow <- writeValues(bflow, do.call(rbind,d$value$value[7,]), bs$row[b])
			tdrain <- writeValues(tdrain, do.call(rbind,d$value$value[9,]), bs$row[b])
			q_in <- writeValues(q_in, do.call(rbind,d$value$value[10,]), bs$row[b])
			
		} else {
			
			matsm[startind[b]:endind[b],] <- do.call(rbind,d$value$value[1,])
			matro[startind[b]:endind[b],] <- do.call(rbind,d$value$value[2,])
			matpet[startind[b]:endind[b],] <- do.call(rbind,d$value$value[3,])
			mataet[startind[b]:endind[b],] <- do.call(rbind,d$value$value[4,])
			matswe[startind[b]:endind[b],] <- do.call(rbind,d$value$value[5,])
			# matcond[startind[b]:endind[b],] <- do.call(rbind,d$value$value[6,])
			matbflow[startind[b]:endind[b],] <- do.call(rbind,d$value$value[7,])
			mattdrain[startind[b]:endind[b],] <- do.call(rbind,d$value$value[9,])
			matq_in[startind[b]:endind[b],] <- do.call(rbind,d$value$value[10,])
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
		# cond <- writeStop(cond)
		bflow <- writeStop(bflow)
		tdrain <- writeStop(tdrain)
		q_in <- writeStop(q_in)		
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
		# cond<-setValues(cond,matcond)
		# baseflow
		bflow<-setValues(bflow,matbflow)
		# time drainage
		tdrain<-setValues(tdrain,mattdrain)
		q_in<-setValues(q_in,matq_in)
		
	}
	close(pb)
	gc()
	return(list(wn=stack(sm),ro=stack(ro),pet=stack(pet),aet=stack(aet),snow=stack(swe),bflow=stack(bflow),tdrain=stack(tdrain),q_in=stack(q_in)))
		
}
