

spacetime = function( p, DATA, OUT=NULL, overwrite=NULL, DS=NULL, method="inla" ) {
	#\\ localized modelling of space and time data to predict/interpolate upon a grid OUT

  p = spacetime.db( p=p, DS="bigmemory.filenames" )
  
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
  
  # set up the data and problem using bigmemory data objects
  print( paste( "Temporary files are being created at:", p$tmp.datadir ) )
  dir.create( p$rootdir, showWarnings=FALSE, recursive =TRUE)
  
  if (is.null(overwrite)) {
    tmpfiles = spacetime.db( p=p, DS="bigmemory.filelist" )
	  for (tf in tmpfiles) {
  		if (file.exists( tf)) {
  			cat( "Temporary files exist from a previous run found. \n")
  			cat( "Send explicit overwrite=TRUE or overwrite=FALSE to proceed. \n") 
  			stop()
  		}
  	}
	}

	if (method=="spatial.covariance" ) {
		if (is.null(overwrite) || overwrite) {
	    spacetime.db( p=p, DS="bigmemory.inputs.data", B=DATA ) # init input data objects
      spacetime.db( p=p, DS="statistics.bigmemory.initialize" ) # init output data objects
		}
	  # define boundary polygon for data .. 
	  # this trims the prediction/statistics locations to speed things up a little ..
	  spacetime.db( p, DS="boundary.redo" ) # ~ 5 min
      if (0) {
        # to reset results manually .. just a template
        # p = spacetime.db( p=p, DS="bigmemory.filenames" )
        S = bigmemory::attach.big.matrix(p$descriptorfile.S, path=p$tmp.datadir)  # statistical outputs
        hist(S[,1] )
        o = which( S[,1] > xxx ) ; S[o,] = NA
        S[sS$problematic,] = NA
        o = which( S[,1] < yyy ); S[o,] = NA
      }
    sS = spacetime.db( p, DS="statistics.bigmemory.status" )
    sS$n.incomplete / ( sS$n.problematic + sS$n.incomplete + sS$n.complete)
    p = make.list( list( jj=sample( sS$incomplete ) ), Y=p ) # random order helps use all cpus
    parallel.run( spacetime.covariance.spatial, p=p ) # no more GMT dependency! :)
    # spacetime.covariance.spatial( p=p )  # if testing serial process
    # save to file
    print( paste( "Results are being saved to:", p$fn.results.covar ) )
    stats = bigmemory::attach.big.matrix(p$descriptorfile.S, path=p$tmp.datadir)  # statistical outputs
    stats = as.data.frame( stats[] )
    save(stats, file=p$fn.results.covar, compress=TRUE )
    print( paste( "Temporary files are being deleted at:", p$tmp.datadir, "tmp" ) )
    spacetime.db( p=p, DS="bigmemory.cleanup" )
    return( p )
  }


  # ------------------------

  
  if (  method =="inla.interpolations" ) {
  	if (is.null(overwrite) || overwrite) {
        spacetime.db( p=p, DS="bigmemory.inputs.data", B=DATA )
        spacetime.db( p=p, DS="bigmemory.inputs.prediction", B=OUT) # covas on prediction locations
        spacetime.db( p=p, DS="statistics.bigmemory.initialize" )
        spacetime.db( p=p, DS="predictions.bigmemory.initialize" )
  	}
    cat( "Warning this will take a very *long* time! (weeks) /n")
    sS = spacetime.db( p, DS="statistics.bigmemory.status" )
    sS$n.incomplete / ( sS$n.problematic + sS$n.incomplete + sS$n.complete)
    p = make.list( list( jj=sample( sS$incomplete ) ), Y=p ) # random order helps use all cpus
    parallel.run( spacetime.interpolate.inla.local, p=p ) # no more GMT dependency! :)
    # spacetime.interpolate.inla.local( p=p, debugrun=TRUE )  # if testing serial process
    if (0) {
      # for checking status of outputs **during** parallel runs: they access bigmemory temporary files
#        bathymetry.figures( DS="statistics", p=p )
#        bathymetry.figures( DS="predictions", p=p )
#        bathymetry.figures( DS="predictions.error", p=p )
      p = spacetime.db( p=p, DS="bigmemory.filenames" )
      S = bigmemory::attach.big.matrix(p$descriptorfile.S, path=p$tmp.datadir)  # statistical outputs
      hist(S[,1] )
      o = which( S[,1] > 600 ); S[o,] = NA
      S[sS$problematic,] = NA
      o = which( S[,1] < 10 );  S[o,] = NA
    }

    # save to file
    pp = bigmemory::attach.big.matrix(p$descriptorfile.P, path=p$tmp.datadir)  # predictions
    preds = pp[]
    ppl = bigmemory::attach.big.matrix(p$descriptorfile.Ploc, path=p$tmp.datadir)
    predloc = ppl[]
    preds = as.data.frame( cbind ( predloc, preds ) )
    names(preds) = c( "plon", "plat", "ndata", "mean", "sdev" )
    save( preds, file=p$fn.P, compress=TRUE )
    rm(preds, predloc)

    # this also rescales results to the full domain
    #\\ statistics are stored at a different resolution than the final grid
    #\\   this fast interpolates the solutions to the final grid
    S = bigmemory::attach.big.matrix(p$descriptorfile.S, path=p$tmp.datadir)  # statistical outputs
    ss = as.data.frame( S[] )
    statnames  = c( "range", "range.sd", "spatial.sd", "observation.sd" )
    statnames0 = c( "range", "range.sd", "spatial.var", "observation.var"  )
    datalink   = c( I(log), I(log), I(log), I(log))   # a log-link seems appropriate for these data
    revlink   = c( I(exp), I(exp), I(exp), I(exp))   # a log-link seems appropriate for these data
    names(ss) = statnames0
    ssl = bigmemory::attach.big.matrix(p$descriptorfile.Sloc, path=p$tmp.datadir)  # statistical output locations
    sslocs = as.data.frame(ssl[]) # copy
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

    # clean up bigmemory files
    spacetime.db( p=p, DS="bigmemory.cleanup" )
    return(p)
  }


}



