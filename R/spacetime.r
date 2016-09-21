

spacetime = function( p, DATA, OUT=NULL, overwrite=NULL, DS=NULL, method=NULL ) {
  #\\ localized modelling of space and time data to predict/interpolate upon a grid OUT

  p = spacetime.db( p=p, DS="filenames" )
  
	if (!is.null(DS)) {
   	#\\ data extraction layer for user objects created by spacetime
    if ( DS=="spatial.covariance") {
      stats = NULL
      if (file.exists( p$fn.results.covar) ) load( p$fn.results.covar )
      return(stats)
    }
    if ( DS =="inla.predictions" ) {
      preds = NULL
      if (file.exists( p$fn.P ) ) load( p$fn.P )
      return( preds )
    }
    if ( DS =="inla.statistics" ) {
      stats = NULL
      if (file.exists( p$fn.S) ) load( p$fn.S )
      return( stats )
    }
	}
  
  # set up the data and problem using data objects
  print( paste( "Temporary files are being created at:", p$tmp.datadir ) )
  dir.create( p$rootdir, showWarnings=FALSE, recursive =TRUE)
  
  if (is.null(overwrite)) {
    tmpfiles = unlist( p$ptr)
	  for (tf in tmpfiles) {
  		if (file.exists( tf)) {
  			cat( "Temporary files exist from a previous run found at: \n")
        cat( tf )
        cat( "\n")
  			cat( "Send an explicit overwrite=TRUE (i.e., restart) or overwrite=FALSE (i.e., continue) option to proceed. \n") 
  			return(NULL)
  		}
  	}
	}

  # init input data
  if (is.null(overwrite) || overwrite) {
    p = spacetime.db( p=p, DS="data.initialize", B=DATA ) # p is updated with pointers to ff data
  }


  if (0) {
    # DEBUG:: for checking status of outputs **during** parallel runs: they access temporary files
    S = p$ff$S  # statistical outputs
    hist(S[,1] )
    o = which( S[,1] > 600 ); S[o,] = NA
    S[sS$problematic,] = NA
    o = which( S[,1] < 10 );  S[o,] = NA
    close(S)
  }


  if (method=="space.time.seasonal" ) { 
    # 2D space and time,  no covariates, eg, temperature
    p$statsvars =  c("varZ", "varSpatial", "varObs", "range", "phi", "kappa", ....) 
    
    if (is.null(overwrite) || overwrite) {
      p = spacetime.db( p=p, DS="predictions.initialize.xyts" )
      p = spacetime.db( p=p, DS="statistics.initialize" ) # init output data objects
    }
    # define boundary polygon for data .. zz a little ..
    if (p$spacetime.stats.boundary.redo) spacetime.db( p, DS="boundary.redo" ) # ~ 5 min
    o = spacetime.db( p, DS="statistics.status" )
    p = make.list( list(jj=sample(  o$incomplete )) , Y=p ) # random order helps use all cpus
    parallel.run( spacetime.interpolate.xyts, p=p ) # no more GMT dependency! :)
    # spacetime.interpolate.xyt( p=p )  # if testing serial process
    # save to file
    print( paste( "Results are being saved to:", p$fn.results.covar ) )
    
    # statistical outputs
    stats = as.data.frame( p$ff$S )
    save(stats, file=p$fn.results.covar, compress=TRUE )
    close(p$ff$S)
    
    print( paste( "Temporary files are being deleted at:", p$tmp.datadir, "tmp" ) )
    spacetime.db( p=p, DS="cleanup" )
    return( p )
  }


	# no time .. pure space, no covariates and no prediction
	if (method=="spatial.covariance" ) { 
    # not used here but passed onto "statistics.initialize" to determine size of output stats matrix
    p$statsvars =  c("varZ", "varSpatial", "varObs", "range", "phi", "kappa") 
		
    if (is.null(overwrite) || overwrite) {
      p = spacetime.db( p=p, DS="statistics.initialize" ) # init output data objects
		}
	  # define boundary polygon for data .. zz a little ..
	  if (p$spacetime.stats.boundary.redo) spacetime.db( p, DS="boundary.redo" ) # ~ 5 min
    o = spacetime.db( p, DS="statistics.status" )
    p = make.list( list(jj=sample(  o$incomplete )) , Y=p ) # random order helps use all cpus
    parallel.run( spacetime.covariance.spatial, p=p ) # no more GMT dependency! :)
    # spacetime.covariance.spatial( p=p )  # if testing serial process
    # save to file
    print( paste( "Results are being saved to:", p$fn.results.covar ) )
    
    # statistical outputs
    stats = as.data.frame( p$ff$S[] )
    save(stats, file=p$fn.results.covar, compress=TRUE )
    close( p$ff$S )

    print( paste( "Temporary files are being deleted at:", p$tmp.datadir, "tmp" ) )
    spacetime.db( p=p, DS="cleanup" )
    return( p )
  }


  # ------------------------

  
  if (  method =="inla.interpolations" ) {
    
    p$statsvars = c("varSpatial", "varObs", "range", "range.sd" )
    nstatvars = length( p$statsvars )

  	# no time .. pure spatial effects and covariates .. 
  	if (is.null(overwrite) || overwrite) {
      p = spacetime.db( p=p, DS="statistics.initialize" ) # init output data objects
      p = spacetime.db( p=p, DS="predictions.initialize.xy", B=OUT )
  	}
    cat( "Warning this will take a very *long* time! (weeks) /n")
    o = spacetime.db( p, DS="statistics.status" )
    p = make.list( list(jj=sample(  o$incomplete )) , Y=p ) # random order helps use all cpus
    parallel.run( spacetime.interpolate.inla.local, p=p ) # no more GMT dependency! :)
    # spacetime.interpolate.inla.local( p=p, debugrun=TRUE )  # if testing serial process
    
    # save to file
    # predictions
    preds = as.data.frame( cbind ( p$ff$Ploc[], p$ff$P[] ) )
    names(preds) = c( "plon", "plat", "ndata", "mean", "sdev" )
    save( preds, file=p$fn.P, compress=TRUE )
    
    rm(preds)
    close(p$ff$Ploc)
    close(p$ff$P)
    # this also rescales results to the full domain
    #\\ statistics are stored at a different resolution than the final grid
    #\\   this fast interpolates the solutions to the final grid
    
    # statistical outputs
    datalink   = c( I(log), I(log), I(log), I(log))   # a log-link seems appropriate for these data
    revlink   = c( I(exp), I(exp), I(exp), I(exp))   # a log-link seems appropriate for these data

    ss = as.data.frame( cbind( p$ff$Sloc[], p$ff$S[] ) )
    names(ss) = c( p$variables$LOCS, p$statsvars )
    
    close(p$ff$S)
    close(p$ff$Sloc)

    locsout = expand.grid( p$plons, p$plats ) # final output grid
    attr( locsout , "out.attrs") = NULL
    names( locsout ) = p$variables$LOCS
    stats = matrix( NA, ncol=nstatvars, nrow=nrow( locsout) )  # output data
    colnames(stats)=p$statsvars
    for ( i in 1:nstatvars ) {
      data = list( x=p$sbbox$plons, y=p$sbbox$plats, z=datalink[[i]](ss[,i]) )
      res = spacetime.interpolate.xy( interp.method="kernel.density", 
        data=ss, locsout=locsout, nr=length(p$plons), nc=length( p$plats),  
        theta=p$dist.mwin, xwidth=p$dist.mwin*10, ywidth=p$dist.mwin*10)
      if (!is.null(res)) stats[i,] = revlink[[i]] (res)
    }
    save( stats,  file=p$fn.S, compress=TRUE )

    plotdata=FALSE ## to debug
    if (plotdata) {
      p$spatial.domain="canada.east"  # force isobaths to work in levelplot
      datarange = log( c( 5, 1200 ))
      dr = seq( datarange[1], datarange[2], length.out=150)
      oc = landmask( db="worldHires", regions=c("Canada", "US"),
                     return.value="not.land", tag="predictions" )  ## resolution of "predictions" which is the final grid size
      toplot = cbind( locsout, z=(stats[,"range"]) )[oc,]
      resol = c(p$dist.mwin,p$dist.mwin)
      levelplot( log(z) ~ plon + plat, toplot, aspect="iso", at=dr, col.regions=color.code( "seis", dr) ,
        contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE), cex=2, resol=resol,
        panel = function(x, y, subscripts, ...) {
          panel.levelplot (x, y, subscripts, aspect="iso", rez=resol, ...)
          cl = landmask( return.value="coast.lonlat",  ylim=c(36,53), xlim=c(-72,-45) )
          cl = lonlat2planar( data.frame( cbind(lon=cl$x, lat=cl$y)), proj.type=p$internal.crs )
          panel.xyplot( cl$plon, cl$plat, col = "black", type="l", lwd=0.8 )
        }
      )
      p$spatial.domain="canada.east.highres"
    }

    # clean up tmp files
    spacetime.db( p=p, DS="cleanup" )
    return(p)
  }

  # -------------------f

  if (  method =="xyts" ) {
     
    p$statsvars = c("varSpatial", "varObs", "range", "range.sd" )
    nstatvars = length( p$statsvars )

    # space, time and covars 
  	if (is.null(overwrite) || overwrite) {
        p = spacetime.db( DS="inputs.prediction", B=OUT) # covas on prediction locations
        p = spacetime.db( DS="statistics.initialize", 
          B=matrix( NA, nrow=p$sbbox$nrow, ncol=p$sbbox$ncol ) )
  	}

  	fmla = as.formula(p$gam.harmonic.formula)  # covariates  

  	# full.model = gam( ). ... -- no random effects

  	# residuals
    o = spacetime.db( p, DS="statistics.status" )
    p = make.list( list(jj=sample(  o$incomplete )) , Y=p ) # random order helps use all cpus
    parallel.run( spacetime.interpolate.xyts, p=p ) # no more GMT dependency! :)
    # spacetime.interpolate.xyts( p=p, debugrun=TRUE )  # if testing serial process

    # save to file
    preds = as.data.frame( cbind ( p$ff$Ploc[], p$ff$P[] ) )
    names(preds) = c( "plon", "plat", "ndata", "mean", "sdev" )
    save( preds, file=p$fn.P, compress=TRUE )
    rm(preds)
    close(p$ff$P)
    close(p$ff$Ploc)
    
    # this also rescales results to the full domain
    #\\ statistics are stored at a different resolution than the final grid
    #\\   this fast interpolates the solutions to the final grid
    datalink   = c( I(log), I(log), I(log), I(log))   # a log-link seems appropriate for these data
    revlink   = c( I(exp), I(exp), I(exp), I(exp))   # a log-link seems appropriate for these data
    
    ss = as.data.frame( cbind( p$ff$Sloc[], p$ff$S[] ) )
    names(ss) = c(p$variables$LOCS, p$statvars)
    close(p$ff$S)
    close(p$ff$Sloc)

    ss$spatial.sd = sqrt( ss$spatial.var )
    ss$observation.sd = sqrt( ss$observation.var )
    ss$spatial.var = NULL
    ss$observation.var = NULL
    locsout = expand.grid( p$plons, p$plats ) # final output grid
    attr( locsout , "out.attrs") = NULL
    names( locsout ) = p$variables$LOCS
    stats = matrix( NA, ncol=nstatvars, nrow=nrow( locsout) )  # output data
    colnames(stats)= p$statsvars
    for ( i in 1:nstatvars ) {
      data = list( x=p$sbbox$plons, y=p$sbbox$plats, z=datalink[[i]](ss[,i]) )
      res = spacetime.interpolate.xy( interp.method="kernel.density", 
        data=ss, locsout=locsout, nr=length(p$plons), nc=length( p$plats),  
        theta=p$dist.mwin, xwidth=p$dist.mwin*10, ywidth=p$dist.mwin*10)
      if (!is.null(res)) stats[i,] = revlink[[i]] (res)
    }
    save( stats,  file=p$fn.S, compress=TRUE )

    plotdata=FALSE ## to debug
    if (plotdata) {
      p$spatial.domain="canada.east"  # force isobaths to work in levelplot
      datarange = log( c( 5, 1200 ))
      dr = seq( datarange[1], datarange[2], length.out=150)
      oc = landmask( db="worldHires", regions=c("Canada", "US"),
                     return.value="not.land", tag="predictions" )  ## resolution of "predictions" which is the final grid size
      toplot = cbind( locsout, z=(stats[,"range"]) )[oc,]
      resol = c(p$dist.mwin,p$dist.mwin)
      levelplot( log(z) ~ plon + plat, toplot, aspect="iso", at=dr, col.regions=color.code( "seis", dr) ,
        contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE), cex=2, resol=resol,
        panel = function(x, y, subscripts, ...) {
          panel.levelplot (x, y, subscripts, aspect="iso", rez=resol, ...)
          cl = landmask( return.value="coast.lonlat",  ylim=c(36,53), xlim=c(-72,-45) )
          cl = lonlat2planar( data.frame( cbind(lon=cl$x, lat=cl$y)), proj.type=p$internal.crs )
          panel.xyplot( cl$plon, cl$plat, col = "black", type="l", lwd=0.8 )
        }
      )
      p$spatial.domain="canada.east.highres"
    }

    # clean up tempfiles
    spacetime.db( p=p, DS="cleanup" )
    return(p)
  }

}



