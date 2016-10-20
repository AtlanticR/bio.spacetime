
spacetime = function( p, DATA, family=gaussian, method="simple", overwrite=NULL, storage.backend="bigmemory.ram" ) {
  #\\ localized modelling of space and time data to predict/interpolate upon a grid OUT
  #\\ overwrite = FALSE restarts from a saved state

  if(0) {
     p = bio.temperature::temperature.parameters( current.year=2016 )
     p = bio.spacetime::spacetime_db( p=p, DS="load.parameters" ) 
     RLibrary( p$libs )
     o = spacetime_db( p, DS="statistics.status" )
     p = make.list( list( locs=sample( o$incomplete )) , Y=p ) 
     spacetime_interpolate (p=p ) 
  }

  p$libs = unique( c( p$libs, "gstat", "sp", "rgdal", "parallel", "mgcv", "fields" ) ) 
  
  p$storage.backend = storage.backend
  if (any( grepl ("ff", p$storage.backend)))         p$libs = c( p$libs, "ff", "ffbase" )
  if (any( grepl ("bigmemory", p$storage.backend)))  p$libs = c( p$libs, "bigmemory" )

  RLibrary( p$libs )
  
  p$spacetime_method = method
  p$spacetime_family = family
  if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )  # default

  p = spacetime_db( p=p, DS="filenames" )
  p$ptr = list() # location for data pointers
  
  # set up the data and problem using data objects
  if (is.null(overwrite)) {
    tmpfiles = unlist( p$diskcache)
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

  if ( is.null(overwrite) || overwrite ) {
  
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

    # number of time slices
    if (!exists("nt", p)) {
      p$nt = 1  
      if (exists( "ny", p)) p$nt = p$nt * p$ny  # annual time slices
      if (exists( "nw", p)) p$nt = p$nt * p$nw  # sub-annual time slices
    } 

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
      p$statsvars = c("varSpatial", "varObs", "range", "range.sd" )# not used .. just for posterity
    }

    message( "Initializing temporary storage of data and outputs (will take a bit longer on NFS clusters) ... ")
    message( "These are large files (4GB each), esp. prediction grids (5 min .. faster if on fileserver), so be patient. ")
    spacetime_db( p=p, DS="cleanup" )
    p = spacetime_db( p=p, DS="statistics.initialize" ) # init output data objects
    p = spacetime_db( p=p, DS="data.initialize", B=DATA$input ) # p is updated with pointers to data
    p = spacetime_db( p=p, DS="predictions.initialize", B=DATA$output )
    
    # p = spacetime_db( p=p, DS="model.covariates.redo", B=DATA$input ) # first pass to model covars only

    rm(DATA); gc()

    message( "Defining boundary polygon for data .. this reduces the number of points to analyse") 
    message( "but takes a few minutes to set up ...")
    spacetime_db( p, DS="boundary.redo" ) # ~ 5 min on nfs
    p = spacetime_db( p, DS="data.filter" )

    spacetime_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to data objects
    message( "Finished. Moving onto analysis... ")

  } else {

    p = spacetime_db( p=p, DS="load.parameters" )  # ie. restart with saved parameters
    RLibrary( p$libs )

  }

  # -------------------------------------
  # 1. localized space-time interpolation
  
  o = spacetime_db( p, DS="statistics.status" )
  p = make.list( list( locs=sample( o$incomplete )) , Y=p ) # random order helps use all cpus
  parallel.run( spacetime_interpolate, p=p ) 

  # 2. fast/simple spatial interpolation for anything not resolved by the local analysis
  p = make.list( list( tiyr_index=p$ncP ), Y=p ) 
  parallel.run( spacetime_interpolate_xy_simple_multiple, p=p ) 
  
  spacetime_db( p, DS="spacetime.predictions.redo" ) # save to disk
  spacetime_db( p, DS="stats.to.prediction.grid.redo")
  
  message ("Finished! \n")
  resp = readline( "To delete temporary files, type <Yes>:  ")
  if (resp=="Yes") {
    spacetime_db( p=p, DS="cleanup" )
  } else {
    message( "Leaving temporary files alone in case you need to examine them or restart a process.")
    message( "You can delete them by running: spacetime_db( p=p, DS='cleanup' ), once you are done.") 
  }
  
  return( p )
  
}



