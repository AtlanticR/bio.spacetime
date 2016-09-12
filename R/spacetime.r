

spacetime = function( p, DATA, OUT=NULL, overwrite=NULL, DS=NULL, method="inla" ) {
  #\\ localized modelling of space and time data to predict/interpolate upon a grid OUT

  p = spacetime.db( p=p, DS="filenames" )
  
	if (!is.null(DS)) {
   	#\\ data extraction layer for user objects created by spacetime
    if ( DS=="spatial.covariance")  {
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
  
  # set up the data and problem using hdf5 data objects
  print( paste( "Temporary files are being created at:", p$tmp.datadir ) )
  dir.create( p$rootdir, showWarnings=FALSE, recursive =TRUE)
  
  if (is.null(overwrite)) {
    tmpfiles = spacetime.db( p=p, DS="filelist" )
	  for (tf in tmpfiles) {
  		if (file.exists( tf)) {
  			cat( "Temporary files exist from a previous run found. \n")
  			cat( "Send explicit overwrite=TRUE or overwrite=FALSE to proceed. \n") 
  			stop()
  		}
  	}
	}

	# no time .. pure space, no covariates and no prediction
	if (method=="spatial.covariance" ) { 
		if (is.null(overwrite) || overwrite) {
			spacetime.db( p=p, DS="data.initialize", B=DATA )
      spacetime.db( p=p, DS="statistics.initialize" ) # init output data objects
		}
	  # define boundary polygon for data .. zz a little ..
	  spacetime.db( p, DS="boundary.redo" ) # ~ 5 min
      if (0) {
        # to reset results manually .. just a template
        # p = spacetime.db( p=p, DS="filenames" )
        S = h5file( p$ptr$S )["S"]  # statistical outputs
        hist(S[,1] )
        o = which( S[,1] > xxx ) ; S[o,] = NA_real_
        S[sS$problematic,] = NA
        o = which( S[,1] < yyy ); S[o,] = NA_real_
        h5close(S)
      }
    p = make.list( list( jj=sample( spacetime.db( p, DS="statistics.status" )["incomplete"] ) ), Y=p ) # random order helps use all cpus
    parallel.run( spacetime.covariance.spatial, p=p ) # no more GMT dependency! :)
    # spacetime.covariance.spatial( p=p )  # if testing serial process
    # save to file
    print( paste( "Results are being saved to:", p$fn.results.covar ) )
    stats = h5file( p$ptr$S)["S"][]  # statistical outputs
    stats = as.data.frame( stats )
    save(stats, file=p$fn.results.covar, compress=TRUE )
    h5close(stats)
    print( paste( "Temporary files are being deleted at:", p$tmp.datadir, "tmp" ) )
    spacetime.db( p=p, DS="cleanup" )
    return( p )
  }


  # ------------------------

  
  if (  method =="inla.interpolations" ) {
  	# no time .. pure spatial effects and covariates .. 
  	if (is.null(overwrite) || overwrite) {
      spacetime.db( p=p, DS="data.initialize", B=DATA )
      spacetime.db( p=p, DS="statistics.initialize" ) # init output data objects
      spacetime.db( p=p, DS="predictions.initialize", B=OUT )
  	}
    cat( "Warning this will take a very *long* time! (weeks) /n")
    
    p = make.list( list( jj=sample( spacetime.db( p, DS="statistics.status" )["incomplete"] ) ), Y=p ) # random order helps 
    parallel.run( spacetime.interpolate.inla.local, p=p ) # no more GMT dependency! :)
    # spacetime.interpolate.inla.local( p=p, debugrun=TRUE )  # if testing serial process
    
    if (0) {
      # for checking status of outputs **during** parallel runs: they access hdf5 temporary files
#        bathymetry.figures( DS="statistics", p=p )
#        bathymetry.figures( DS="predictions", p=p )
#        bathymetry.figures( DS="predictions.error", p=p )
      p = spacetime.db( p=p, DS="filenames" )
      S = h5file( p$ptr$S)["S"]  # statistical outputs
      hist(S[,1] )
      o = which( S[,1] > 600 ); S[o,] = NA_real_
      S[sS$problematic,] = NA
      o = which( S[,1] < 10 );  S[o,] = NA_real_
      h5close(S)
    }

    # save to file
    pp = h5file( p$ptr$P)["P"]  # predictions
    preds = pp[]
    h5close(pp)
    
    ppl = h5file( p$ptr$Ploc)["Ploc"]
    predloc = ppl[]
    h5close(ppl)

    preds = as.data.frame( cbind ( predloc, preds ) )
    names(preds) = c( "plon", "plat", "ndata", "mean", "sdev" )
    save( preds, file=p$fn.P, compress=TRUE )
    
    rm(preds, predloc)

    # this also rescales results to the full domain
    #\\ statistics are stored at a different resolution than the final grid
    #\\   this fast interpolates the solutions to the final grid
    S = h5file( p$ptr$S)["S"]  # statistical outputs
    ss = as.data.frame( S[] )
    h5close(S)

    statnames  = c( "range", "range.sd", "spatial.sd", "observation.sd" )
    statnames0 = c( "range", "range.sd", "spatial.var", "observation.var"  )
    datalink   = c( I(log), I(log), I(log), I(log))   # a log-link seems appropriate for these data
    revlink   = c( I(exp), I(exp), I(exp), I(exp))   # a log-link seems appropriate for these data
    names(ss) = statnames0
    
    ssl = h5file( p$ptr$Sloc)["Sloc"]  # statistical output locations
    sslocs = as.data.frame(ssl[]) # copy
    h5close(ssl)
    
    names(sslocs) = p$variables$LOCS
    ss = cbind( sslocs, ss )
    rm (S)
    ss$spatial.sd = sqrt( ss$spatial.var )
    ss$observation.sd = sqrt( ss$observation.var )
    ss$spatial.var = NULL
    ss$observation.var = NULL
    locsout = expand.grid( p$plons, p$plats ) # final output grid
    attr( locsout , "out.attrs") = NULL
    names( locsout ) = p$variables$LOCS
    stats = matrix( NA, ncol=length(statnames), nrow=nrow( locsout) )  # output data
    colnames(stats)=statnames
    for ( i in 1:length(statnames) ) {
      data = list( x=p$sbbox$plons, y=p$sbbox$plats, z=datalink[[i]](ss[,i]) )
      res = spacetime.reshape( data=ss, locsout=locsout, nr=length(p$plons), nc=length( p$plats),  
        interp.method="kernel.density", theta=p$dist.mwin, nsd=10)
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

    # clean up hdf5 files
    spacetime.db( p=p, DS="cleanup" )
    return(p)
  }

  # -------------------f

  if (  method =="gam.harmonic" ) {
  	# space, time and covars 
  	if (is.null(ovehdf5| overwrite) {
        spacetime.db( DS="inputs.data", B=DATA )
        spacetime.db( DS="inputs.prediction", B=OUT) # covas on prediction locations
        spacetime.db( DS="statistics.initialize", B=matrix( NA_real_, nrow=p$sbbox$nrow, ncol=p$sbbox$ncol ) )
  	}

  	fmla = as.formula(p$gam.harmonic.formula)  # covariates  

  	full.model = gam( ). ... -- no random effects

  	residuals
  	
    p = make.list( list( jj=sample( spacetime.db( p, DS="statistics.status" )["incomplete"] ) ), Y=p ) # random order helps 
    parallel.run( spacetime.interpolate.gam.harmonic, p=p ) # no more GMT dependency! :)
    # spacetime.interpolate.gam.harmonic( p=p, debugrun=TRUE )  # if testing serial process
    

    if (0) {
      # for checking status of outputs **during** parallel runs: they access hdf5 temporary files
#        bathymetry.figures( DS="statistics", p=p )
#        bathymetry.figures( DS="predictions", p=p )
#        bathymetry.figures( DS="predictions.error", p=p )
      p = spacetime.db( p=p, DS="filenames" )
      S = h5file( p$ptr$S)["S"]  # statistical outputs
      hist(S[,1] )
      o = which( S[,1] > 600 ); S[o,] = NA
      S[sS$problematic,] = NA
      o = which( S[,1] < 10 );  S[o,] = NA
      h5close(S)
    }

    # save to file
    pp = h5file( p$ptr$P)["P"]  # predictions
    preds = pp[]
    h5close(pp)
    
    ppl = h5file( p$ptr$Ploc)["Ploc"]
    predloc = ppl[]
    h5close(ppl)
    
    preds = as.data.frame( cbind ( predloc, preds ) )
    names(preds) = c( "plon", "plat", "ndata", "mean", "sdev" )
    save( preds, file=p$fn.P, compress=TRUE )
    rm(preds, predloc)

    # this also rescales results to the full domain
    #\\ statistics are stored at a different resolution than the final grid
    #\\   this fast interpolates the solutions to the final grid
    S = h5file( p$ptr$S)["S"]  # statistical outputs
    ss = as.data.frame( S[] )
    h5close(S)
    
    statnames  = c( "range", "range.sd", "spatial.sd", "observation.sd" )
    statnames0 = c( "range", "range.sd", "spatial.var", "observation.var"  )
    datalink   = c( I(log), I(log), I(log), I(log))   # a log-link seems appropriate for these data
    revlink   = c( I(exp), I(exp), I(exp), I(exp))   # a log-link seems appropriate for these data
    names(ss) = statnames0
    
    ssl = h5file( p$ptr$Sloc)["Sloc"]  # statistical output locations
    sslocs = as.data.frame(ssl[]) # copy
    h5close(ssl)
    names(sslocs) = p$variables$LOCS
    
    ss = cbind( sslocs, ss )
    rm (S)
    ss$spatial.sd = sqrt( ss$spatial.var )
    ss$observation.sd = sqrt( ss$observation.var )
    ss$spatial.var = NULL
    ss$observation.var = NULL
    locsout = expand.grid( p$plons, p$plats ) # final output grid
    attr( locsout , "out.attrs") = NULL
    names( locsout ) = p$variables$LOCS
    stats = matrix( NA, ncol=length(statnames), nrow=nrow( locsout) )  # output data
    colnames(stats)=statnames
    for ( i in 1:length(statnames) ) {
      data = list( x=p$sbbox$plons, y=p$sbbox$plats, z=datalink[[i]](ss[,i]) )
      res = spacetime.reshape( data=ss, locsout=locsout, nr=length(p$plons), nc=length( p$plats),  
        interp.method="kernel.density", theta=p$dist.mwin, nsd=10)
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

    # clean up hdf5 files
    spacetime.db( p=p, DS="cleanup" )
    return(p)
  }



}



