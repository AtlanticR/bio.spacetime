

spacetime = function( p, DATA, overwrite=NULL, method=NULL ) {
  #\\ localized modelling of space and time data to predict/interpolate upon a grid OUT

  p$libs = unique( c( p$libs, "gstat", "sp", "rgdal", "parallel", "mgcv", "ff", "ffbase", "fields" ) )

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

  # permit passing a function rather than data directly .. less RAM usage
  if (class(DATA)=="character") assign("DATA", eval(parse(text=DATA) ) )


  # -------------------------------------

  if (method=="xy" ) { 

    if ( p$spacetime_engine=="inla" ) {
      # kept for historical reasons ..
      p$statsvars = c("varSpatial", "varObs", "range", "range.sd" )
      nstatvars = length( p$statsvars )
      if (!exists("modelformula", p) ) {
        p$modelformula = " -1 + intercept + f( spatial.field, model=SPDE ) " # SPDE is the spatial object to be created by inla 
        p$modelformula = as.formula( paste( p$variables$Y, "~", p$modelformula ) )
      }
      if (is.null(overwrite) || overwrite) {
        p = spacetime_db( p=p, DS="statistics.initialize" ) # init output data objects: requires p$statvars
        p = spacetime_db( p=p, DS="data.initialize", B=DATA$input ) # p is updated with pointers to ff data
        p = spacetime_db( p=p, DS="predictions.initialize.xy", B=DATA$output )
        rm(DATA);gc()
        spacetime_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to ff objectss
        message( "Finshed. Moving onto analysis... ")
      } else {
        p = spacetime_db( p=p, DS="load.parameters" )
      }
      message( "Warning this will take a very *long* time! (weeks) /n")
      o = spacetime_db( p, DS="statistics.status" )
      p = make.list( list( locs=sample(  o$incomplete )) , Y=p ) # random order helps use all cpus
      parallel.run( spacetime_interpolate_xy_local_inla, p=p ) # no more GMT dependency! :)
      # save to file
      preds = as.data.frame( cbind ( p$ff$Ploc[], p$ff$P[] ) )
      names(preds) = c( "plon", "plat", "ndata", "mean", "sdev" )
      save( preds, file=p$fn.P, compress=TRUE )
      rm(preds)
      close(p$ff$Ploc)
      close(p$ff$P)
      # this also rescales results to the full domain
      datalink   = c( I(log), I(log), I(log), I(log))   # a log-link seems appropriate for these data
      revlink   = c( I(exp), I(exp), I(exp), I(exp))   # a log-link seems appropriate for these data
      ss = as.data.frame( cbind( p$ff$Sloc[], p$ff$S[] ) )
      names(ss) = c( p$variables$LOCS, p$statsvars )
      close(p$ff$S)
      close(p$ff$Sloc)
      locsout = expand.grid( p$plons, p$plats ) # final output grid
      attr( locsout , "out.attrs") = NULL
      names( locsout ) = p$variables$LOCS
      stats = matrix( NaN, ncol=nstatvars, nrow=nrow( locsout) )  # output data
      colnames(stats)=p$statsvars
      for ( i in 1:nstatvars ) {
        data = list( x=p$sbbox$plons, y=p$sbbox$plats, z=datalink[[i]](ss[,i]) )
        res = spacetime_interpolate_xy_singlepass( interp.method="kernel.density", 
          data=ss, locsout=locsout, nr=length(p$plons), nc=length( p$plats),  
          theta=p$dist.mwin, xwidth=p$dist.mwin*10, ywidth=p$dist.mwin*10)
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
      spacetime_db( p=p, DS="cleanup" )
      return(p)
    }

    # ------------
    
    if ( p$spacetime_engine=="gam" ) {
      p$statsvars =  c("varZ", "varSpatial", "varObs", "range", "phi", "kappa") 
      if (exists("modelformula", p) ) {
        #message( "The formula is equivalent to Y ~ 1 + spatial.covariance" ) 
      }
      if (is.null(overwrite) || overwrite) {
        p = spacetime_db( p=p, DS="data.initialize", B=DATA$input ) # p is updated with pointers to ff data
        p = spacetime_db( p=p, DS="statistics.initialize" ) # init output data objects: requires p$statvars
        spacetime_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to ff  
      } else {
        p = spacetime_db( p=p, DS="load.parameters" )
      }
      # define boundary polygon for data .. zz a little ..
      if (p$spacetime.stats.boundary.redo) spacetime_db( p, DS="boundary.redo" ) # ~ 5 min
      o = spacetime_db( p, DS="statistics.status" )
      p = make.list( list( locs=sample( o$incomplete )) , Y=p ) # random order helps use all cpus
      parallel.run( spacetime_interpolate_xy, p=p ) # no more GMT dependency! :)

      #  save predictions to disk

      # statistical outputs
      stats = as.data.frame( cbind( p$ff$Sloc, p$ff$S ) )
      save(stats, file=p$stats, compress=TRUE )
      close( p$ff$S )

      print( paste( "Temporary files are being deleted at:", p$tmp.datadir, "tmp" ) )
      spacetime_db( p=p, DS="cleanup" )
      return( p )
    }
  }



  # -------------------------------------

  if (method=="xyt" ) { 
 
  }


  # -------------------------------------

  if (method=="xyts" ) { 
    # 2D space and time, eg, temperature
    # spacetime_engine specific methods are inside of sapcetime_interpolate_xyts 
    #.. this is req here to identifiy the size of the statistics output dim
    if ( p$spacetime_engine %in% 
      c( "harmonics.1", "harmonics.2", "harmonics.3", "harmonics.1.depth",
         "seasonal.basic", "seasonal.smoothed", "annual"  ) ) {
      p$statsvars =  c("varZ" )  ## to add to this we need to add to the function creating stats in  spacetime_timeseries
    }

    if (is.null(overwrite) || overwrite) {
      message( "Initializing temporary storage of data and outputs ... ")
      message( "These are large files so be patient. ")
      p = spacetime_db( p=p, DS="statistics.initialize" ) # init output data objects
      p = spacetime_db( p=p, DS="data.initialize", B=DATA$input ) # p is updated with pointers to ff data
      p = spacetime_db( p=p, DS="predictions.initialize.xyts", B=DATA$output )
      rm(DATA); gc()
      spacetime_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to ff objectss
      message( "Finshed. Moving onto analysis... ")
    } else {
      p = spacetime_db( p=p, DS="load.parameters" )
    }
    if (p$spacetime.stats.boundary.redo) {
      message( "Defining boundary polygon for data .. this reduces the number of points to analyse") 
      message( "but takes a few (~5) minutes to complete ...")
      spacetime_db( p, DS="boundary.redo" ) # ~ 5 min
    }
    if (!exists("modelformula", p) ) {
      # these are simple, generic defaults .. 
      # for more complex models (.i.e, with covariates) the formula should be passed directly 
      p$modelformula = switch( p$spacetime_engine ,
        seasonal.basic = ' s(yr) + s(dyear, bs="cc") ', 
        seasonal.smoothed = ' s(yr, dyear) + s(yr) + s(dyear, bs="cc")  ', 
        harmonics.1 = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w)  ', 
        harmonics.2 = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) ' , 
        harmonics.3 = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) + s(yr, cos.w3) + s(yr, sin.w3)  + s(cos.w3) + s( sin.w3 ) ',
        harmonics.1.depth = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) +s(z)  ', 
        annual = ' s(yr) '
      )
      p$modelformula = as.formula( paste( p$variables$Y, "~", p$modelformula ) )
      message( "Verify that the modelformula is/should be:" )
      message( p$modelformula )
    }
 
    # 1. localized space-time interpolation
    o = spacetime_db( p, DS="statistics.status" )
    p = make.list( list( locs=sample( o$incomplete )) , Y=p ) # random order helps use all cpus
    parallel.run( spacetime_interpolation_xyts, p=p ) 

    # 2. fast/simple spatial interpolation for anything not yet resolved
    # .. but first, pre-compute a few things 
    Ploc = p$ff$Ploc
    p$Mat2Ploc = cbind( (Ploc[,1]-p$plons[1])/p$pres + 1, (Ploc[,2]-p$plats[1])/p$pres + 1) # row, col indices in matrix form
    close(Ploc)
    p$wght = setup.image.smooth(nrow=p$nplons, ncol=p$nplats, dx=p$pres, dy=p$pres,
      theta=p$theta, xwidth=p$nsd*p$theta, ywidth=p$nsd*p$theta )
    p = make.list( list( tiyr_index=1:(p$nw*p$yr)), Y=p ) # random order helps use all cpus
    parallel.run( spacetime_interpolate_xy_simple_multiple, p=p ) 

    # 3. save to disk
    ffP   = p$ff$P  # predictions
    ffPsd = p$ff$Psd  # predictions
    for ( r in 1:length(p$tyears) ) {
      y = p$tyears[r]
      fn1 = file.path( p$savedir, paste("spacetime.interpolation",  y, "rdata", sep="." ) )
      fn2 = file.path( p$savedir, paste("spacetime.interpolation.sd",  y, "rdata", sep="." ) )
      col.ranges = (r-1) * p$nw + (1:p$nw) 
      P = ffP  [,col.ranges]
      V = ffPsd[,col.ranges]
      save( P, file=fn1, compress=T )
      save( V, file=fn2, compress=T )
      print ( paste("Year:", y)  )
    }
    close(ffP)
    close(ffPsd)

    # statistical outputs
    stats = as.data.frame( cbind( p$ff$Sloc, p$ff$S ) )
    ## warp this onto prediction grid
    save(stats, file=p$fn.stats, compress=TRUE )
    close(p$ff$S)

    ## Here we parse predictions and save in a more useful format

    print( paste( "Temporary files are being deleted at:", p$tmp.datadir, "tmp" ) )
    spacetime_db( p=p, DS="cleanup" )
    return( p )
  }
    
}



