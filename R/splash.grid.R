#' Simple process-led algorithms for simulating habitats (SPLASH v.2.0)
#' 
#' R/C++ implementation of the SPLASH v.2.0 algorithm (Davis et al., 2017; Sandoval et al., in prep.).
#' 
#' @param   sw_in Incoming shortwave solar radiation (W m-2), Raster* object of monthly or daily averages with z time dimension.
#' @param   tc Air temperature (°C), same dimensions as sw_in
#' @param   pn Precipitation (mm), same dimensions as sw_in
#' @param   elev Elevation (m.a.s.l)
#' @param   soil Raster* object with the layers organized as sand(perc),clay(perc),organic matter(perc),coarse-fragments-fraction(perc), bulk density(g cm-3) and depth(m)
#' @param   outdir (optional) directory path where the results will be saved, the working directory by default
#' @param   tmpdir (optional) directory path where the temporary files will be saved
#' @param   sim.control (optional) list including options to control the output:
#' - monthly_out = TRUE by default; daily if FALSE.
#' - inmem = FALSE by default; writes all the results to the disk by chunks to save RAM,
#'       sacrificing speed.
#' @return a list of rasterBricks objects with z time dimension, all of them saved to outdir as netcdf files:
#' \itemize{
#'         \item \eqn{wn}: Soil water content (mm) within the first 2 m of depth.
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
#' \dontrun{data(SA_cru)
#' splash.grid(
#'   sw_in = SA_cru$sw_in,
#'   tc = SA_cru$tc,
#'   pn = SA_cru$pn,
#'   elev = SA_cru$elev,
#'   soil = SA_cru$soil,
#'   outdir = getwd(),
#'   tmpdir = dirname(raster::rasterTmpFile()),
#'   sim.control = list(
#'     monthly_out = FALSE,
#'     inmem = TRUE
#'   )
#' )}

splash.grid<-function(sw_in, tc, pn, elev, soil, outdir=getwd(),tmpdir=dirname(rasterTmpFile()),sim.control=list(monthly_out=TRUE,inmem=FALSE)){
	###########################################################################
	# 00. Check if parallel computation is required by the user and if the dimensions of the raster objects match
	###########################################################################
	on.exit(endCluster())
	clcheck<-try(getCluster(), silent=TRUE)
	if(class(clcheck)[1]=="try-error"){
		# If no cluster is initialized, assume only one core will do the calculations, beginCluster(1) saved me the time of coding serial versions of the functions
		beginCluster(1,'SOCK')
		message('Only using one core, use first beginCluster(ncores) if you want to run splash in parallel!!')
		
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
			# xmat <- raster::getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
			# ymat <- raster::getValues(y, row=bs$row[i], nrows=bs$nrows[i] )
			xncells<-cellFromRow(x,bs$row[i]:(bs$row[i]+ bs$nrows[i]-1))
			#xmat<-raster::getValues(x,bs$row[i], bs$nrows[i])
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
	if (ncell(elev)>1e7|(res(elev)[1]>0.1 & res(elev)[1]<1)){
		cat("computing contributing areas and terrain features, it might take a while...","\n")
		start.time<-Sys.time()
		Au<-upslope_areav2(elev, type,tmpdir)
		# 1.2.2 get latitude from the high res dem
		lat<-getlatitude(elev, filename='lat.grd')
		system("gdaldem slope -s 111120 -compute_edges rawdem.tif slope_deg.tif")
		system("gdaldem aspect -zero_for_flat -compute_edges rawdem.tif aspect_deg.tif")
		terraines<-raster::stack(list(slope='slope_deg.tif',aspect='aspect_deg.tif'))
		# get resolution in m2
		resolution<-area(elev,filename='area.grd',overwrite=T)
		resolution<-calc(resolution,function(x){sqrt(x)*1000})
		gc()
		
		end.time<-Sys.time()
		timetaken<-end.time-start.time
		cat(paste0("... it took"),timetaken,"\n")
		
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
		##fix slopes NA as flats
		terraines[is.na(terraines) & !is.na(elev)]<-0
		cat("done!","\n")
		gc()
	}
		
	
	#lat<-writeRaster(lat,filename="lat.grd", overwrite=TRUE)
	
	# 1.3  calculate slope and aspect
	
	
	###########################################################################
	# 02. Get time info from data
	###########################################################################
	### ztime: the time index of the inputs
	### ztime_out: the time index of the outputs
	################################################
	##BCE dates https://github.com/RMHoek/NOAAearthquakeAnalysis/blob/master/R/as_BC_date.R
	################################################
	as_BC_date <- function(year, month = 1, day = 1){
		if(year < 0) year<-(-year)
		Y <- as.character(year)
		M <- as.character(month)
		D <- as.character(day)
		fwdY <- paste(Y, "1", "1", sep = "/")
		fwdYMD <- paste(Y, M, D, sep = "/")
		AD0 <- lubridate::as_date("0000/1/1") ##merry xmas!
		n_AD0 <- as.numeric(AD0)
		n_fwdY <- as.numeric(lubridate::as_date(fwdY))
		n_MD <- as.numeric(lubridate::as_date(fwdYMD)) -
		as.numeric(lubridate::as_date(fwdY))
		n_BC <- n_AD0 - (n_fwdY - n_AD0) + n_MD
		if(n_MD==0) n_BC <- n_BC + 1
		BC_date <- lubridate::as_date(n_BC)
		return(BC_date)
	}
	
	
	ztime<-getZ(pn)
	
	if(is.character(ztime[1])){
		timechar<-do.call(rbind,strsplit(ztime,'-',fixed = T))
		y<-as.numeric(unique(timechar[,1]))
		ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
		ztime1<-as_BC_date(y[1], month = 1, day = 1)
		ztime1<-as_BC_date(timechar[1,1], month = timechar[1,2], day = timechar[1,3])
		ztime2<-as_BC_date(timechar[2,1], month = timechar[2,2], day = timechar[2,3])
		time.freq<-ztime1-ztime2
		ztime.months<-seq(as_BC_date(y[1], month = 1, day = 1),as_BC_date(y[length(y)], month = 12, day = 31), by="month")
		ztime.days<-seq(as_BC_date(y[1], month = 1, day = 1),as_BC_date(y[length(y)], month = 12, day = 31), by="days")
		zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31))
		
	}else if (is(ztime[1],'Date')){
		y<-as.numeric(unique(format(ztime,'%Y')))
		ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
		
		ztime.months<-seq(as.Date(paste(y[1],1,sep="-"),format="%Y-%j"),as.Date(paste(y[length(y)],ny[length(y)],sep="-"),format="%Y-%j"), by="month")
		ztime.days<-seq(as.Date(paste(y[1],1,sep="-"),format="%Y-%j"),as.Date(paste(y[length(y)],ny[length(y)],sep="-"),format="%Y-%j"), by="day")
		zunit=paste("days","since",paste0(y[1]-1,"-",12,"-",31))
		time.freq<-ztime[2]-ztime[1]
	}
	###########################################################################
	# 03. create empty rasters
	###########################################################################
	### define the number of layers
	
	if(sim.control$monthly_out){
		time_text<-'monthly'
		nl=length(ztime.months)
		ztime_out <- ztime.months
	}else{
		time_text<-'daily'
		nl=length(ztime.days)
		ztime_out <- ztime.days
	}
	
	
	setwd(outdir)
		
	if(!sim.control$inmem){
		
		
		# actual soil moisture
		sm<-brick(elev, values=FALSE, nl=nl,filename=paste0(outdir,"/",'SPLASH_v2.0_',"swc",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="wm", varunit="mm",	longname='Soil water content', xname="lon", yname="lat", zname="time")
		sm<-setZ(sm,ztime_out)
		
		ro<-brick(elev, values=FALSE, nl=nl,filename=paste0(outdir,"/",'SPLASH_v2.0_',"ro",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="ro", varunit="mm",	longname='surface runoff', xname="lon", yname="lat", zname="time")
		ro<-setZ(ro,ztime_out)
		
		pet<-brick(elev, values=FALSE, nl=nl,filename=paste0(outdir,"/",'SPLASH_v2.0_',"pet",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="pet", varunit="mm",	longname='potential evapotranspiration', xname="lon", yname="lat", zname="time")
		pet<-setZ(pet,ztime_out)
		
		aet<-brick(elev, values=FALSE, nl=nl,filename=paste0(outdir,"/",'SPLASH_v2.0_',"aet",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="aet", varunit="mm",	longname='actual evapotranspiration', xname="lon", yname="lat", zname="time")
		aet<-setZ(aet,ztime_out)
		
		
		swe<-brick(elev, values=FALSE, nl=nl,filename=paste0(outdir,"/",'SPLASH_v2.0_',"swe",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="swe", varunit="mm",	longname='snow water equivalent', xname="lon", yname="lat", zname="time")
		swe<-setZ(swe,ztime_out)
		
		cond<-brick(elev, values=FALSE, nl=nl,filename=paste0(outdir,"/",'SPLASH_v2.0_',"cond",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="cond", varunit="mm",	longname='condensation', xname="lon", yname="lat", zname="time")
		cond<-setZ(cond,ztime_out)
		
		bflow<-brick(elev, values=FALSE, nl=nl,filename=paste0(outdir,"/",'SPLASH_v2.0_',"bflow",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="bflow", varunit="mm",	longname='drainage', xname="lon", yname="lat", zname="time")
		bflow<-setZ(bflow,ztime_out)
		
		netrad<-brick(elev, values=FALSE, nl=nl,filename=paste0(outdir,"/",'SPLASH_v2.0_',"netrad",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="netrad", varunit="MJ",	longname='daytime net radiation', xname="lon", yname="lat", zname="time")
		netrad<-setZ(netrad,ztime_out)
		
		sm_lim<-brick(elev, values=FALSE, nl=nl,filename=paste0(outdir,"/",'SPLASH_v2.0_',"sm_lim",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="sm_lim", varunit="mm",	longname='fraction of available water content', xname="lon", yname="lat", zname="time")
		sm_lim<-setZ(sm_lim,ztime_out)
		
		
		
	}else{
		# actual soil moisture
		sm<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=nl)
		extent(sm)<-extent(elev)
		# runoff
		ro<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=nl)
		extent(ro)<-extent(elev)
		# Potential evapotranspiration
		pet<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=nl)
		extent(pet)<-extent(elev)
		# Actual evapotranspiration
		aet<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=nl)
		extent(aet)<-extent(elev)
		# Snow water equivalent
		swe<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=nl)
		extent(swe)<-extent(elev)
		# Condensation
		cond<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=nl)
		extent(cond)<-extent(elev)
		# baseflow
		bflow<-brick(nrows=nrow(elev), ncols=ncol(elev), crs=crs(elev), nl=nl)
		extent(bflow)<-extent(elev)
		netrad<-bflow
		sm_lim<-bflow
		#move inputs to the memory
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
	nodes <- length(cl)
	########################### set the size of the blocks
	message('Using cluster with ', nodes, ' nodes')
	if(nodes>1){
		bs <- blockSize(sw_in, minblocks=nodes*2)
	}else{
		bs <- blockSize(sw_in, minblocks=nodes*10)
	}
	########################### export the vaiables to the nodes 
	parallel::clusterExport(cl, c("sw_in","tc","pn","elev","lat","terraines",'soil','resolution','Au','ztime','bs','splash.point','sim.control'),envir=environment()) 
	pb <- pbCreate(bs$n)
	pb <- txtProgressBar(min=1,max = max(bs$n,2), style = 3)
	#cat("computing...","\n")
	###############################################################################################
	# 02. create the functions to send to the workers, split the data in chunks
	###############################################################################################	
	clFun <- function(i) {
		swrow<- split(raster::getValues(sw_in, bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		tcrow<-split(raster::getValues(tc,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		pnrow<-split(raster::getValues(pn,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		elevrow<-raster::getValues(elev,bs$row[i], bs$nrows[i])
		sloprow<-raster::getValues(terraines[[1]],bs$row[i], bs$nrows[i])
		asprow<-raster::getValues(terraines[[2]],bs$row[i], bs$nrows[i])
		latrow<-raster::getValues(lat,bs$row[i], bs$nrows[i])
		soilrow<-split(raster::getValues(soil,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		Aurow<-split(raster::getValues(Au,bs$row[i], bs$nrows[i]),1:(ncol(elev)*bs$nrows[i]))
		resrow<-raster::getValues(resolution,bs$row[i], bs$nrows[i])
		# do calculations
		
		
		splash_block<-mapply(splash.point,
			sw_in=swrow,	# shortwave radiation W/m2
			tc=tcrow,		# air temperature C
			pn= pnrow,		# precipitation mm
			lat=latrow,		# latitude deg
			elev=elevrow,		# elevation masl
			slop=sloprow,	# slope deg
			asp=asprow,		# aspect deg
			soil_data=soilrow, 		
			Au=Aurow,
			resolution=resrow,
			MoreArgs=list(ts_out=FALSE,monthly_out=sim.control$monthly_out,time_index=ztime,verbose=FALSE),
			SIMPLIFY =TRUE
		)
		
		
		return(splash_block)
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
	if(!sim.control$inmem){
		
		
		sm<-writeStart(sm,filename=paste0(outdir,"/",'SPLASH_v2.0_',"sm",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="wm", varunit="mm",	longname='Soil water content', xname="lon", yname="lat", zname="time")
		
		ro<-writeStart(ro,filename=paste0(outdir,"/",'SPLASH_v2.0_',"ro",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="ro", varunit="mm",	longname='surface runoff', xname="lon", yname="lat", zname="time")
		
		pet<-writeStart(pet,filename=paste0(outdir,"/",'SPLASH_v2.0_',"pet",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="pet", varunit="mm",	longname='potential evapotranspiration', xname="lon", yname="lat", zname="time")
		
		aet<-writeStart(aet,filename=paste0(outdir,"/",'SPLASH_v2.0_',"aet",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="aet", varunit="mm",	longname='actual evapotranspiration', xname="lon", yname="lat", zname="time")
		
		swe<-writeStart(swe,filename=paste0(outdir,"/",'SPLASH_v2.0_',"swe",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="swe", varunit="mm",	longname='snow water equivalent', xname="lon", yname="lat", zname="time")
		
		cond<-writeStart(cond,filename=paste0(outdir,"/",'SPLASH_v2.0_',"cond",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="cond", varunit="mm",	longname='condensation', xname="lon", yname="lat", zname="time")
		
		bflow<-writeStart(bflow,filename=paste0(outdir,"/",'SPLASH_v2.0_',"bflow",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="bflow", varunit="mm",	longname='drainage', xname="lon", yname="lat", zname="time")
		
		netrad<-writeStart(netrad,filename=paste0(outdir,"/",'SPLASH_v2.0_',"netrad",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="netrad", varunit="MJ",	longname='daytime net radiation', xname="lon", yname="lat", zname="time")
		
		sm_lim<-writeStart(sm_lim,filename=paste0(outdir,"/",'SPLASH_v2.0_',"sm_lim",'_',time_text,'_',y[1],"-",y[length(y)],".nc"),format="CDF",overwrite=TRUE,varname="sm_lim", varunit="mm",	longname='fraction of available water content', xname="lon", yname="lat", zname="time")
		
		
		
		
		
	}else {
		matsm <- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matro<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matpet <- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		mataet <- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matswe<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matcond <- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matbflow<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matnetrad<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		matsm_lim<- matrix(ncol=nlayers(sm), nrow=ncell(sm))
		endind<-cumsum(bs$nrows*sm@ncols)
		startind<-c(1,endind+1)    
	}
	###############################################################################################
	# 05. receive results from the nodes
	###############################################################################################	
	for (i in 1:bs$n) {
		
		d <- parallel::recvOneData(cl)
		# error?
		if (! d$value$success) {
			stop('cluster error:',"\n",d$value$value)
		}
		# which block is this?
		b <- d$value$tag
		# cat('received block: ',b,'\n'); flush.console();
		if (!sim.control$inmem) {
			sm <- writeValues(sm,do.call(rbind,d$value$value[1,]), bs$row[b])
			ro <- writeValues(ro, do.call(rbind,d$value$value[2,]), bs$row[b])
			pet <- writeValues(pet, do.call(rbind,d$value$value[3,]), bs$row[b])
			aet <- writeValues(aet, do.call(rbind,d$value$value[4,]), bs$row[b])
			swe <- writeValues(swe, do.call(rbind,d$value$value[5,]), bs$row[b])
			cond <- writeValues(cond, do.call(rbind,d$value$value[6,]), bs$row[b])
			bflow <- writeValues(bflow, do.call(rbind,d$value$value[7,]), bs$row[b])
			netrad <- writeValues(netrad, do.call(rbind,d$value$value[8,]), bs$row[b])
			sm_lim <- writeValues(sm_lim, do.call(rbind,d$value$value[9,]), bs$row[b])
			
			
		} else {
			
			matsm[startind[b]:endind[b],] <- do.call(rbind,d$value$value[1,])
			matro[startind[b]:endind[b],] <- do.call(rbind,d$value$value[2,])
			matpet[startind[b]:endind[b],] <- do.call(rbind,d$value$value[3,])
			mataet[startind[b]:endind[b],] <- do.call(rbind,d$value$value[4,])
			matswe[startind[b]:endind[b],] <- do.call(rbind,d$value$value[5,])
			matcond[startind[b]:endind[b],] <- do.call(rbind,d$value$value[6,])
			matbflow[startind[b]:endind[b],] <- do.call(rbind,d$value$value[7,])
			matnetrad[startind[b]:endind[b],] <- do.call(rbind,d$value$value[8,])
			matsm_lim[startind[b]:endind[b],] <- do.call(rbind,d$value$value[9,])
			}
		
		# need to send more data?
		ni <- nodes + i
		if (ni <= bs$n) {
			parallel:::sendCall(cl[[d$node]], clFun, ni, tag=ni)
		}
		setTxtProgressBar(pb,i)
	}
	close(pb)
	###############################################################################################
	# 06. close connection with the files, or assign valueas to the raster objects
	###############################################################################################
	
	if (!sim.control$inmem) {
		sm <- writeStop(sm)
		ro <- writeStop(ro)
		pet <- writeStop(pet)
		aet <- writeStop(aet)
		swe <- writeStop(swe)
		cond <- writeStop(cond)
		bflow <- writeStop(bflow)
		netrad <- writeStop(netrad)
		sm_lim <- writeStop(sm_lim)		
		
	} else {
		# soil water content
		sm<-setValues(sm,matsm)
		sm<-setZ(sm,ztime_out)
		# runoff
		ro<-setValues(ro,matro)
		ro<-setZ(ro,ztime_out)
		# Potential evapotranspiration
		pet<-setValues(pet,matpet)
		pet<-setZ(pet,ztime_out)
		# Actual evapotranspiration
		aet<-setValues(aet,mataet)
		aet<-setZ(aet,ztime_out)
		# Snow water equivalent
		swe<-setValues(swe,matswe)
		swe<-setZ(swe,ztime_out)
		# Condensation
		cond<-setValues(cond,matcond)
		cond<-setZ(cond,ztime_out)
		# baseflow
		bflow<-setValues(bflow,matbflow)
		bflow<-setZ(bflow,ztime_out)
		# net radiation
		netrad<-setValues(netrad,matnetrad)
		netrad<-setZ(netrad,ztime_out)
		##fraction FC
		sm_lim<-setValues(sm_lim,matsm_lim)
		sm_lim<-setZ(sm_lim,ztime_out)
		
	}
	
	gc()
	return(list(wn=sm,ro=ro,pet=pet,aet=aet,snow=swe,bflow=bflow,cond=cond,netrad=netrad,sm_lim=sm_lim))
		
		
}

