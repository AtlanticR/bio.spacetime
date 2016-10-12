
spacetime = function( p, DATA, overwrite=NULL) {
  #\\ localized modelling of space and time data to predict/interpolate upon a grid OUT

  p$libs = unique( c( p$libs, "gstat", "sp", "rgdal", "parallel", "mgcv", "bigmemory", "fields" ) )

  if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )  # default
  
  p$tmp.datadir = file.path( p$project.root, "tmp" )
  message( paste( "Temporary files are being created at:", p$tmp.datadir ) )
  if( !file.exists(p$tmp.datadir)) dir.create( p$tmp.datadir, recursive=TRUE, showWarnings=FALSE )

  p$savedir = file.path(p$project.root, "spacetime", p$spatial.domain )
  
  message( paste( "Final outputs will be palced at:", p$savedir ) )
  if( !file.exists(p$savedir)) dir.create( p$savedir, recursive=TRUE, showWarnings=FALSE )
    
  p = spacetime_db( p=p, DS="filenames" )
  
  # set up the data and problem using data objects
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

  if (!exists("spacetime_engine_modelformula", p) ) {
    # these are simple, generic defaults .. 
    # for more complex models (.i.e, with covariates) the formula should be passed directly 
    p$spacetime_engine_modelformula = switch( p$spacetime_engine ,
      seasonal.basic = ' s(yr) + s(dyear, bs="cc") ', 
      seasonal.smoothed = ' s(yr, dyear) + s(yr) + s(dyear, bs="cc")  ', 
      harmonics.1 = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w)  ', 
      harmonics.2 = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) ' , 
      harmonics.3 = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) + s(yr, cos.w3) + s(yr, sin.w3)  + s(cos.w3) + s( sin.w3 ) ',
      harmonics.1.depth = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) +s(z)  ', 
      inla = ' -1 + intercept + f( spatial.field, model=SPDE ) ', # not used
      annual = ' s(yr) ',
      '+1'  # default, ie. Y~ + 1 , just to catch error
    )
    p$spacetime_engine_modelformula = as.formula( paste( p$variables$Y, "~", p$spacetime_engine_modelformula ) )
    message( "Verify that the spacetime_engine_modelformula is/should be:" )
    message( p$spacetime_engine_modelformula )
  }


  # permit passing a function rather than data directly .. less RAM usage
  if (class(DATA)=="character") assign("DATA", eval(parse(text=DATA) ) )

  # require knowledge of size of stats output before create S, which varies with a given type of analysis


  if ( p$spacetime_engine %in% c( "harmonics.1", "harmonics.2", "harmonics.3", "harmonics.1.depth",
         "seasonal.basic", "seasonal.smoothed", "annual", "gam"  ) ) {
    p$statsvars = c( "sdTotal", "rsquared", "ndata" )
  }
  
  if (exists("spacetime_variogram_engine", p) ) {
    p$statsvars = c( p$statsvars, "sdSpatial", "sdObs", "range", "phi", "nu")
  }

  if ( exists("TIME", p$variables) ){
    p$statsvars = c( p$statsvars, "ar_timerange", "ar_1" )
  }

  if ( p$spacetime_engine=="inla") {
    # not used .. just for posterity
    p$statsvars = c("varSpatial", "varObs", "range", "range.sd" )
  }

  if (is.null(overwrite) || overwrite) {
    message( "Initializing temporary storage of data and outputs (will take a bit longer on NFS clusters) ... ")
    message( "These are large files so be patient. ")
    spacetime_db( p=p, DS="cleanup" )
    p = spacetime_db( p=p, DS="statistics.initialize" ) # init output data objects
    p = spacetime_db( p=p, DS="data.initialize", B=DATA$input ) # p is updated with pointers to data
    p = spacetime_db( p=p, DS="model.covariates.redo", B=DATA$input ) # first pass to model covars
    p = spacetime_db( p=p, DS="predictions.initialize", B=DATA$output )
    spacetime_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to data objects
    message( "Finshed. Moving onto analysis... ")
    rm(DATA); gc()
  } else {
    p = spacetime_db( p=p, DS="load.parameters" )
  }

  if (p$spacetime.stats.boundary.redo) {
    message( "Defining boundary polygon for data .. this reduces the number of points to analyse") 
    message( "but takes a few minutes to set up ...")
    spacetime_db( p, DS="boundary.redo" ) # ~ 5 min on nfs
  }

  # -------------------------------------
  # 1. localized space-time interpolation
  o = spacetime_db( p, DS="statistics.status" )
  p = make.list( list( locs=sample( o$incomplete )) , Y=p ) # random order helps use all cpus
  parallel.run( spacetime_interpolate, p=p ) 

  # 2. fast/simple spatial interpolation for anything not yet resolved
  # .. but first, pre-compute a few things 
  Ploc = attach.big.matrix( p$ptr$Ploc )
  p$Mat2Ploc = cbind( (Ploc[,1]-p$plons[1])/p$pres + 1, (Ploc[,2]-p$plats[1])/p$pres + 1) # row, col indices in matrix form
  p$wght = setup.image.smooth( nrow=p$nplons, ncol=p$nplats, dx=p$pres, dy=p$res, 
    theta=p$theta, xwidth=p$nsd*p$theta, ywidth=p$nsd*p$theta )

  # a little more interpolation
  if (exists("TIME", p)) {
    if (exists("nw", p)) {
      p = make.list( list( tiyr_index=1:(p$nw*p$yr)), Y=p ) 
    } else {
      p = make.list( list( tiyr_index=1:p$yr), Y=p ) 
    }
    parallel.run( spacetime_interpolate_xy_simple_multiple, p=p ) 
  } else {
    p = make.list( list( tiyr_index=1), Y=p ) # random order helps use all cpus
    spacetime_interpolate_xy_simple_multiple( p=p )
  }

  spacetime_db( p, DS="spacetime.predictions.redo" ) # save to disk
  spacetime_db( p, DS="stats_to_prediction_grid.redo")

  print( paste( "Temporary files are being deleted at:", p$tmp.datadir, "tmp" ) )
  
  pause( "stop here and check results before cleanuP ...")
  #spacetime_db( p=p, DS="cleanup" )


  return( p )


  ## -- finished -- rest kept for historical reasons
 
  if ( 0 ) {
    # kept for historical reasons .. inla methods
    nstatvars = length( p$statsvars )
    message( "Warning this will take a very *long* time! (weeks) /n")
    o = spacetime_db( p, DS="statistics.status" )
    p = make.list( list( locs=sample(  o$incomplete )) , Y=p ) # random order helps use all cpus
    parallel.run( spacetime_interpolate_xy_local_inla, p=p ) # no more GMT dependency! :)
    # save to file
    P = attach.big.matrix( p$ptr$P )
    Ploc = attach.big.matrix( p$ptr$Ploc )
    preds = as.data.frame( cbind ( Ploc[], P[] ) )
    names(preds) = c( "plon", "plat", "ndata", "mean", "sdev" )
    save( preds, file=p$fn.P, compress=TRUE )
    rm(preds)
    # this also rescales results to the full domain
    datalink   = c( I(log), I(log), I(log), I(log))   # a log-link seems appropriate for these data
    revlink   = c( I(exp), I(exp), I(exp), I(exp))   # a log-link seems appropriate for these data
    S = attach.big.matrix( p$ptr$S )
    Sloc = attach.big.matrix( p$ptr$Sloc )
    ss = as.data.frame( cbind( Sloc[], S[] ) )
    names(ss) = c( p$variables$LOCS, p$statsvars )
    locsout = expand.grid( p$plons, p$plats ) # final output grid
    attr( locsout , "out.attrs") = NULL
    names( locsout ) = p$variables$LOCS
    stats = matrix( NaN, ncol=nstatvars, nrow=nrow( locsout) )  # output data
    colnames(stats)=p$statsvars
    for ( i in 1:nstatvars ) {
      data = list( x=p$sbbox$plons, y=p$sbbox$plats, z=datalink[[i]](ss[,i]) )
      res = spacetime_interpolate_xy_singlepass( interp.method="kernel.density", 
        data=ss, locsout=locsout, nr=length(p$plons), nc=length( p$plats),  
        theta=p$spacetime.prediction.dist.min, xwidth=p$spacetime.prediction.dist.min*10, ywidth=p$spacetime.prediction.dist.min*10)
      if (!is.null(res)) stats[i,] = revlink[[i]] (res)
    }
    save( stats,  file=p$fn.S, compress=TRUE )

    if (0) {
      p$spatial.domain="canada.east"  # force isobaths to work in levelplot
      datarange = log( c( 5, 1200 ))
      dr = seq( datarange[1], datarange[2], length.out=150)
      oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="predictions" )  
      ## resolution of "predictions" which is the final grid size
      toplot = cbind( locsout, z=(stats[,"range"]) )[oc,]
      resol = c(p$spacetime.prediction.dist.min,p$spacetime.prediction.dist.min)
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
    spacetime_db( p=p, DS="cleanup" )
    return(p)
  }
  
}



