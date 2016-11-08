
spacetime = function( p, DATA, family=gaussian(), overwrite=NULL, storage.backend="bigmemory.ram", boundary=TRUE, do.secondstage=TRUE ) {
  #\\ localized modelling of space and time data to predict/interpolate upon a grid OUT
  #\\ overwrite = FALSE restarts from a saved state
  #\\ speed ratings: bigmemory.ram (1), ff (2), bigmemory.filebacked (3)

  p$stloc = file.path( p$project.root, "tmp" )
  # message( paste( "Temporary files are being created at:", p$stloc ) )
  if( !file.exists(p$stloc)) dir.create( p$stloc, recursive=TRUE, showWarnings=FALSE )

  p$savedir = file.path(p$project.root, "spacetime", p$spatial.domain )
  
  # message( paste( "Final outputs will be palced at:", p$savedir ) )
  if( !file.exists(p$savedir)) dir.create( p$savedir, recursive=TRUE, showWarnings=FALSE )
  
  if(0) {
     p = bio.temperature::temperature.parameters( current.year=2016 )
     family=gaussian()
     overwrite=NULL
     storage.backend="bigmemory.ram"
     DATA='hydro.db( p=p, DS="spacetime.input" )'

     p = bio.temperature::temperature.parameters( current.year=2016 )
     p$stloc = file.path( p$project.root, "tmp" )
     p = bio.spacetime::spacetime_db( p=p, DS="load.parameters" ) 
     RLibrary( p$libs )
     o = spacetime_db( p, DS="statistics.status" )
     p = make.list( list( locs=sample( o$todo )) , Y=p ) 
     spacetime_interpolate (p=p ) 
  }

  p$libs = unique( c( p$libs, "gstat", "sp", "rgdal", "parallel", "mgcv", "fields" ) ) 
  
  p$storage.backend = storage.backend
  if (any( grepl ("ff", p$storage.backend)))         p$libs = c( p$libs, "ff", "ffbase" )
  if (any( grepl ("bigmemory", p$storage.backend)))  p$libs = c( p$libs, "bigmemory" )

  RLibrary( p$libs )
  
  
  # family handling copied from glm
  if (!exists( "spacetime_family", p)) {
    if (is.character(family)) 
      stop( "Please send family as a function") 
    if (is.function(family)) 
      family = family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
    p$spacetime_family = family
  }

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
  
  # these are some possible models:
  #   if (!exists("spacetime_engine_modelformula", p) ) {
    # these are simple, generic defaults .. 
    # for more complex models (.i.e, with covariates) the formula should be passed directly 
    #   p$spacetime_engine_modelformula = switch( p$spacetime_engine ,
    #     seasonal.basic = ' s(yr) + s(dyear, bs="cc") ', 
    #     seasonal.smoothed = ' s(yr, dyear) + s(yr) + s(dyear, bs="cc")  ', 
    #     seasonal.smoothed.depth.lonlat = ' s(yr, dyear) + s(yr, k=3) + s(dyear, bs="cc") +s(z) +s(plon) +s(plat) + s(plon, plat, by=yr), s(plon, plat, k=10, by=dyear ) ', 
    #     seasonal.smoothed.depth.lonlat.complex = ' s(yr, dyear, bs="ts") + s(yr, k=3, bs="ts") + s(dyear, bs="cc") +s(z, bs="ts") +s(plon, bs="ts") +s(plat, bs="ts") + s(plon, plat, by=tiyr, k=10, bs="ts" ) ', 
    #     harmonics.1 = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w)  ', 
    #     harmonics.2 = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) ' , 
    #     harmonics.3 = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) + s(yr, cos.w3) + s(yr, sin.w3)  + s(cos.w3) + s( sin.w3 ) ',
    #     harmonics.1.depth = ' s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) +s(z)  ', 
    #     harmonics.1.depth.lonlat = 's(yr, k=5, bs="ts") + s(cos.w, bs="ts") + s(sin.w, bs="ts") +s(z, k=3, bs="ts") +s(plon,k=3, bs="ts") +s(plat, k=3, bs="ts") + s(plon, plat, cos.w, sin.w, yr, k=100, bs="ts") ', 
    #     inla = ' -1 + intercept + f( spatial.field, model=SPDE ) ', # not used
    #     annual = ' s(yr) ',
    #     '+1'  # default, ie. Y~ + 1 , just to catch error
    #   )
    #   p$spacetime_engine_modelformula = as.formula( paste( p$variables$Y, "~", p$spacetime_engine_modelformula ) )
    #   message( "Verify that the spacetime_engine_modelformula is/should be:" )
    #   message( p$spacetime_engine_modelformula )
    # }

    p$variables$ALL = all.vars( p$spacetime_engine_modelformula )
    
    # permit passing a function rather than data directly .. less RAM usage
    if (class(DATA)=="character") assign("DATA", eval(parse(text=DATA) ) )
    
    # number of time slices
    if (!exists("nt", p)) {
      p$nt = 1  # default to 1 == no time
      if (exists( "ny", p)) p$nt = p$nt * p$ny  # annual time slices
      if (exists( "nw", p)) p$nt = p$nt * p$nw  # sub-annual time slices
    } 

    # require knowledge of size of stats output before create S, which varies with a given type of analysis
    p$statsvars = c( "sdTotal", "rsquared", "ndata" )
    if (p$spacetime_engine == "habitat") p$statsvars = c( p$statsvars)
    if (p$spacetime_engine=="inla") p$statsvars = c(p$statsvars, "varSpatial", "varObs", "range", "range.sd" )# not used .. just for posterity
    if (exists("spacetime_variogram_engine", p) ) p$statsvars = c( p$statsvars, "sdSpatial", "sdObs", "range", "phi", "nu")
    if (exists("TIME", p$variables) ) p$statsvars = c( p$statsvars, "ar_timerange", "ar_1" )
   

    message( "Initializing temporary storage of data and outputs (will take a bit longer on NFS clusters) ... ")
    message( "These are large files (4 to 6 X 5GB), esp. prediction grids (5 min .. faster if on fileserver), so be patient. ")
    spacetime_db( p=p, DS="cleanup" )
    

    # NOTE:: must not sink this into a deeper funcion as bigmemory RAM seems to losse the pointers if they are not made simultaneously (at the same namespace depth) ..

    # init output data objects
    # statistics storage matrix ( aggregation window, coords ) .. no inputs required
    
      # statistics coordinates
      Sloc = as.matrix( expand.grid( p$sbbox$plons, p$sbbox$plats ))
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Sloc = big.matrix(nrow=nrow(Sloc), ncol=ncol(Sloc), type="double"  )
          p$bm$Sloc[] = Sloc
          p$ptr$Sloc  = bigmemory::describe( p$bm$Sloc  )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Sloc  = p$cache$Sloc
          bigmemory::as.big.matrix( Sloc, type="double", backingfile=basename(p$bm$Sloc), descriptorfile=basename(p$cache$Sloc), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Sloc = ff( Sloc, dim=dim(Sloc), file=p$cache$Sloc, overwrite=TRUE )
        }

      
      S = matrix( NaN, nrow=nrow(Sloc), ncol=length( p$statsvars ) ) # NA forces into logical
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$S = big.matrix(nrow=nrow(Sloc), ncol=length( p$statsvars ), type="double"  )
          p$bm$S[] = S
          p$ptr$S  = bigmemory::describe( p$bm$S )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$S  = p$cache$S
          bigmemory::as.big.matrix( S, type="double", backingfile=basename(p$bm$S), descriptorfile=basename(p$cache$S), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$S = ff( S, dim=dim(S), file=p$cache$S, overwrite=TRUE )
        }

      
      Sflag = matrix( NaN, nrow=nrow(Sloc), ncol=1 ) 
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Sflag = big.matrix(nrow=nrow(Sloc), ncol=1, type="double" )
          p$bm$Sflag[] = NaN
          p$ptr$Sflag  = bigmemory::describe( p$bm$Sflag )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Sflag  = p$cache$Sflag
          bigmemory::as.big.matrix( Sflag, type="double", backingfile=basename(p$bm$Sflag), descriptorfile=basename(p$cache$Sflag), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Sflag = ff( Sflag, dim=dim(Sflag), file=p$cache$Sflag, overwrite=TRUE )
        }

      rm(S, Sflag, Sloc)


      # dependent variable
      Y = as.matrix(DATA$input[, p$variables$Y ])
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Y = big.matrix( nrow=nrow(Y), ncol=1, type="double"  )
          p$bm$Y[] = Y
          p$ptr$Y  = bigmemory::describe( p$bm$Y )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Y  = p$cache$Y
          bigmemory::as.big.matrix( Y, type="double", backingfile=basename(p$bm$Y), descriptorfile=basename(p$cache$Y), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Y = ff( Y, dim=dim(Y), file=p$cache$Y, overwrite=TRUE )
        }
      rm(Y)

      # limits based on quantiles to permit in predictions 
      Y = spacetime_attach( p$storage.backend, p$ptr$Y )
      p$qs = quantile( Y[], probs=p$quantile_bounds, na.rm=TRUE  )


      if (p$spacetime_engine == "habitat") {
        logitY = spacetime_db( p=p, DS="presence.absense" )
          if (p$storage.backend == "bigmemory.ram" ) {
            if (!exists("habitat.threshold.quantile", p)) p$habitat.threshold.quantile = 0.01 
            p$bm$Ylogit = big.matrix( nrow=nrow(logitY), ncol=1, type="double" )
            p$bm$Ylogit[] = logitY
            p$ptr$Ylogit  = bigmemory::describe( p$bm$Ylogit )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ylogit  = p$cache$Ylogit
            bigmemory::as.big.matrix( logitY, type="double", backingfile=basename(p$bm$Ylogit), descriptorfile=basename(p$cache$Ylogit), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ylogit = ff( Ylogit, dim=dim(logitY), file=p$cache$Ylogit, overwrite=TRUE )
          }
        rm(logitY)        
      }

     # data coordinates
      Yloc = as.matrix( DATA$input[, p$variables$LOCS ])
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Yloc = big.matrix( nrow=nrow(Yloc), ncol=ncol(Yloc), type="double" )
          p$bm$Yloc[] = Yloc
          p$ptr$Yloc = bigmemory::describe( p$bm$Yloc )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Yloc  = p$cache$Yloc
          bigmemory::as.big.matrix( Yloc, type="double", backingfile=basename(p$bm$Yloc), descriptorfile=basename(p$cache$Yloc), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Yloc = ff( Yloc, dim=dim(Yloc), file=p$cache$Yloc, overwrite=TRUE )
        }
      rm(Yloc)
        

      # independent variables/ covariate
      if (exists("COV", p$variables)) {
        Ycov = as.matrix(  DATA$input[ , p$variables$COV ] )
          if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Ycov = big.matrix( nrow=nrow(Ycov), ncol=ncol(Ycov), type="double")
            p$bm$Ycov[] = Ycov  
            p$ptr$Ycov  = bigmemory::describe( p$bm$Ycov )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ycov  = p$cache$Ycov
            bigmemory::as.big.matrix( Ycov, type="double", backingfile=basename(p$bm$Ycov), descriptorfile=basename(p$cache$Ycov), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ycov = ff( Ycov, dim=dim(Ycov), file=p$cache$Ycov, overwrite=TRUE )
          }
        rm(Ycov)
      }


      # data times
      if ( exists("TIME", p$variables) ) {
        Ytime = as.matrix(  DATA$input[, p$variables$TIME ] )
          if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Ytime = big.matrix( nrow=nrow(Ytime), ncol=ncol(Ytime), type="double"  )
            p$bm$Ytime[] = Ytime
            p$ptr$Ytime  = bigmemory::describe( p$bm$Ytime )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ytime  = p$cache$Ytime
            bigmemory::as.big.matrix( Ytime, type="double", backingfile=basename(p$bm$Ytime), descriptorfile=basename(p$cache$Ytime), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ytime = ff( Ytime, dim=dim(Ytime), file=p$cache$Ytime, overwrite=TRUE )
          }
        rm(Ytime)
      }
  

      if (exists("COV", p$variables)) {
        if (is.vector(DATA$output$COV) ) {
          Pcov = as.matrix( DATA$output$COV ) 
        } else {
          Pcov = as.matrix( DATA$output$COV[,p$variables$COV ] ) 
        }
        attr( Pcov, "dimnames" ) = NULL
          if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Pcov = big.matrix( nrow=nrow(Pcov), ncol=ncol(Pcov), type="double"  )
            p$bm$Pcov[] = Pcov
            p$ptr$Pcov  = bigmemory::describe( p$bm$Pcov )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Pcov  = p$cache$Pcov
            bigmemory::as.big.matrix( Pcov, type="double", backingfile=basename(p$bm$Pcov), descriptorfile=basename(p$cache$Pcov), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Pcov = ff( Pcov, dim=dim(Pcov), file=p$cache$Pcov, overwrite=TRUE )
          }
        rm(Pcov)
      }

      # prediction times 
      if (exists("TIME", p$variables)) {
        Ptime = as.matrix( DATA$output$TIME )
        attr( Ptime, "dimnames" ) = NULL
          if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Ptime = big.matrix( nrow=nrow(Ptime), ncol=ncol(Ptime), type="double" )
            p$bm$Ptime[] = Ptime
            p$ptr$Ptime  = bigmemory::describe( p$bm$Ptime )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ptime  = p$cache$Ptime
            bigmemory::as.big.matrix( Ptime, type="double", backingfile=basename(p$bm$Ptime), descriptorfile=basename(p$cache$Ptime), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ptime = ff( Ptime, dim=dim(Ptime), file=p$cache$Ptime, overwrite=TRUE )
          }
        rm(Ptime)
      }

      # predictions and associated stats
      P = matrix( NaN, nrow=nrow(DATA$output$LOCS), ncol=p$nt )
        # predictions
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$P = big.matrix( nrow=nrow(P), ncol=ncol(P), type="double" )
          p$bm$P[] = P
          p$ptr$P  = bigmemory::describe( p$bm$P )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$P  = p$cache$P
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$P), descriptorfile=basename(p$cache$P), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$P = ff( P, dim=dim(P), file=p$cache$P, overwrite=TRUE )
        }
      
      # count of prediction estimates
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Pn = big.matrix( nrow=nrow(P), ncol=ncol(P), type="double" )
          p$bm$Pn[] = P
          p$ptr$Pn = bigmemory::describe( p$bm$Pn )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Pn  = p$cache$Pn
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Pn), descriptorfile=basename(p$cache$Pn), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Pn = ff( P, dim=dim(P), file=p$cache$Pn, overwrite=TRUE )
        }

      # sd of prediction estimates
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Psd = big.matrix( nrow=nrow(P), ncol=ncol(P), type="double" )
          p$bm$Psd[] = P
          p$ptr$Psd =bigmemory::describe( p$bm$Psd )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Psd  = p$cache$Psd
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Psd), descriptorfile=basename(p$cache$Psd), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Psd = ff( P, dim=dim(P), file=p$cache$Psd, overwrite=TRUE )
        }

        if (p$spacetime_engine == "habitat") {
          if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Plogit= big.matrix( nrow=nrow(P), ncol=ncol(P) , type="double" )
            p$bm$Plogit[] = P
            p$ptr$Plogit = bigmemory::describe(p$bm$Plogit )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Plogit  = p$cache$Plogit
            bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Plogit), descriptorfile=basename(p$cache$Plogit), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Plogit = ff( P, dim=dim(P), file=p$cache$Plogit, overwrite=TRUE )
          }

          if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Plogitsd= big.matrix( nrow=nrow(P), ncol=ncol(P) , type="double" )
            p$bm$Plogitsd[] = P
            p$ptr$Plogitsd = bigmemory::describe(p$bm$Plogitsd )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Plogitsd  = p$cache$Plogitsd
            bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Plogitsd), descriptorfile=basename(p$cache$Plogitsd), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Plogitsd = ff( P, dim=dim(P), file=p$cache$Plogitsd, overwrite=TRUE )
          }
        }
      rm(P)

      # prediction coordinates
      Ploc = as.matrix( DATA$output$LOCS )
      attr( Ploc, "dimnames" ) = NULL
         if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$Ploc = big.matrix( nrow=nrow(Ploc), ncol=ncol(Ploc), type="double" )
            p$bm$Ploc[] = Ploc
            p$ptr$Ploc  = bigmemory::describe( p$bm$Ploc )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ploc  = p$cache$Ploc
            bigmemory::as.big.matrix( Ploc, type="double", backingfile=basename(p$bm$Ploc), descriptorfile=basename(p$cache$Ploc), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ploc = ff( Ploc, dim=dim(Ploc), file=p$cache$Ploc, overwrite=TRUE )
          }

        # pre-compute a few things for spacetime_interpolate_xy_simple_multiple  
        Mat2Ploc = as.matrix( cbind( 
          (Ploc[,1]-p$plons[1])/p$pres + 1, 
          (Ploc[,2]-p$plats[1])/p$pres + 1) ) # row, col indices in matrix form

        attr( Mat2Ploc, "dimnames" ) = NULL
            if (p$storage.backend == "bigmemory.ram" ) {
              p$bm$Mat2Ploc = big.matrix( nrow=nrow(Mat2Ploc), ncol=ncol(Mat2Ploc), type="double" )
              p$bm$Mat2Ploc[] = Mat2Ploc
              p$ptr$Mat2Ploc  = bigmemory::describe( p$bm$Mat2Ploc )
            }
            if (p$storage.backend == "bigmemory.filebacked" ) {
              p$ptr$Mat2Ploc  = p$cache$Mat2Ploc
              bigmemory::as.big.matrix( Mat2Ploc, type="double", backingfile=basename(p$bm$Mat2Ploc), descriptorfile=basename(p$cache$Mat2Ploc), backingpath=p$stloc )
            }
            if (p$storage.backend == "ff" ) {
              p$ptr$Mat2Ploc = ff( Mat2Ploc, dim=dim(Mat2Ploc), file=p$cache$Mat2Ploc, overwrite=TRUE )
            }
        rm(Mat2Ploc)
      rm(Ploc)

      p$spatial_weights = setup.image.smooth( nrow=p$nplons, ncol=p$nplats, dx=p$pres, dy=p$pres, 
        theta=p$theta, xwidth=p$nsd*p$theta, ywidth=p$nsd*p$theta )

      if (0) {
        # to add gloabl covariate model ?? 
        P0   = matrix( 0, nrow=nrow(DATA$output$LOCS), ncol=p$nt )
          if (p$storage.backend == "bigmemory.ram" ) {
             p$bm$P0 = big.matrix( nrow=nrow(P0), ncol=ncol(P0), type="double" )
             p$bm$P0[] = P0
             p$ptr$P0  = bigmemory::describe( p$bm$P0 )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$P0  = p$cache$P0
            bigmemory::as.big.matrix( P0, type="double", backingfile=basename(p$bm$P0), descriptorfile=basename(p$cache$P0), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$P0 = ff( P0, dim=dim(P0), file=p$cache$P0, overwrite=TRUE )
          }

        # P0sd

          if (p$storage.backend == "bigmemory.ram" ) {
            p$bm$P0sd = big.matrix( nrow=nrow(P0), ncol=ncol(P0), type="double" )
            p$bm$P0sd[] = P0
            p$ptr$P0sd  = bigmemory::describe( p$bm$P0sd )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$P0sd  = p$cache$P0sd
            bigmemory::as.big.matrix( P0, type="double", backingfile=basename(p$bm$P0sd), descriptorfile=basename(p$cache$P0sd), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$P0sd = ff( P0, dim=dim(P0), file=p$cache$P0sd, overwrite=TRUE )
          }
        rm(P0)
      }

    if (exists("model.covariates.globally", p) && p$model.covariates.globally ) {
      p = spacetime_db( p=p, DS="model.covariates.redo", B=DATA$input ) # first pass to model covars only
    }

    rm(DATA); gc()

    if (boundary) {
      message( "Defining boundary polygon for data .. this reduces the number of points to analyse") 
      message( "but takes a few minutes to set up ...")
      spacetime_db( p, DS="boundary.redo" ) # ~ 5 min on nfs
    # last set of filters to reduce problem size
      Sflag = spacetime_attach( p$storage.backend, p$ptr$Sflag )
      bnds = try( spacetime_db( p, DS="boundary" ) )
      if (!is.null(bnds)) {
        if( !("try-error" %in% class(bnds) ) ) {
          to.ignore = which( bnds$inside.polygon == 0 ) # outside boundary
          if (length(to.ignore)>0) Sflag[to.ignore,] = Inf
      }}
      bnds = NULL
    }

      Y = spacetime_attach( p$storage.backend, p$ptr$Y )
      Yloc = spacetime_attach( p$storage.backend, p$ptr$Yloc )

      Yi = 1:length(Y) # index with useable data
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) Yi[bad] = NA

      # data locations
      bad = which( !is.finite( rowSums(Yloc[])))
      if (length(bad)> 0 ) Yi[bad] = NA

    # data locations
      if (exists("COV", p$variables)) {
        Ycov = spacetime_attach( p$storage.backend, p$ptr$Ycov )
        if (length(p$variables$COV)==1) {
          bad = which( !is.finite( Ycov[] ))
        } else {
          bad = which( !is.finite( rowSums(Ycov[])))
        }
        if (length(bad)> 0 ) Yi[bad] = NA
        Yi = na.omit(Yi)
      }

      # data locations
      if (exists("TIME", p$variables)) {
        Ytime = spacetime_attach( p$storage.backend, p$ptr$Ytime )
        bad = which( !is.finite( Ytime[] ))
        if (length(bad)> 0 ) Yi[bad] = NA
        Yi = na.omit(Yi)
      }
      bad = NULL

      Yi = as.matrix(Yi)
        if (p$storage.backend == "bigmemory.ram" ) {
          p$bm$Yi = big.matrix( nrow=nrow(Yi), ncol=ncol(Yi), type="double" )
          p$bm$Yi[] = Yi
          p$ptr$Yi  = bigmemory::describe( p$bm$Yi )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Yi  = p$cache$Yi
          bigmemory::as.big.matrix( Yi, type="double", backingfile=basename(p$bm$Yi), descriptorfile=basename(p$cache$Yi), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Yi = ff( Yi, dim=dim(Yi), file=p$cache$Yi, overwrite=TRUE )
        }
      rm(Yi)


      #---------------------
      # prediction locations and covariates
      Ploc = spacetime_attach( p$storage.backend, p$ptr$Ploc )
      p$rcP = data.frame( cbind( 
        Prow = (Ploc[,1]-p$plons[1])/p$pres + 1,  
        Pcol = (Ploc[,2]-p$plats[1])/p$pres + 1) )
      p$rcP$rc = paste( p$rcP$Prow, p$rcP$Pcol, sep="~")
      p$rcP$Prow = p$rcP$Pcol = NULL
      
      #-----------------
      # row, col indices
      # statistical output locations
      Sloc = spacetime_attach( p$storage.backend, p$ptr$Sloc )
      p$rcS = data.frame( cbind( 
        Srow = (Sloc[,1]-p$plons[1])/p$pres + 1,  
        Scol = (Sloc[,2]-p$plats[1])/p$pres + 1))

    # misc intermediate calcs to be done outside of parallel loops
    p$dist.median = (p$dist.max + p$dist.min ) / 2
    p$upsampling = sort( p$sampling[ which( p$sampling > 1 ) ] )
    p$upsampling = p$upsampling[ which(p$upsampling*p$dist.median <= p$dist.max )]
    p$downsampling = sort( p$sampling[ which( p$sampling < 1) ] , decreasing=TRUE )
    p$downsampling = p$downsampling[ which(p$downsampling*p$dist.median >= p$dist.min )]

    spacetime_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to data objects
    message( "Finished. Moving onto analysis... ")
    gc()

  } else {

    p = spacetime_db( p=p, DS="load.parameters" )  # ie. restart with saved parameters
    RLibrary( p$libs )

  }


  # -------------------------------------
  # localized space-time modelling/interpolation/prediction

  o = spacetime_db( p, DS="statistics.status" )
  p = make.list( list( locs=sample( o$todo )) , Y=p ) # random order helps use all cpus
  p$time.start =  Sys.time()
  parallel.run( spacetime_interpolate, p=p ) 
  p$time.end1 =  Sys.time()
  message( paste( "Time taken:", difftime( p$time.end1, p$time.start ) ) )

  # save solutions to disk before continuuing
  spacetime_db( p, DS="spacetime.predictions.redo" ) # save to disk for use outside spacetime*
  spacetime_db( p, DS="stats.to.prediction.grid.redo") # save to disk for use outside spacetime*
  
  if ( do.secondstage ) {
    # 2. same interpolation method but relax the spatial extent
    o = spacetime_db( p, DS="statistics.reset.problem.locations" )
    if (length(o$todo) > 0) {
      p$spacetime_distance_prediction = p$spacetime_distance_prediction * 2
      p$dist.max = p$dist.max * 2 
      p = make.list( list( locs=sample( o$todo )) , Y=p ) # random order helps use all cpus
      parallel.run( spacetime_interpolate, p=p ) 
      p$time.end2 =  Sys.time()
      message( paste( "Time taken to stage 2:", difftime( p$time.end2, p$time.end1 ) ) )
    }
  }
  # save solutions to disk (again .. overwrite)
    spacetime_db( p, DS="spacetime.predictions.redo" ) # save to disk for use outside spacetime*
    spacetime_db( p, DS="stats.to.prediction.grid.redo") # save to disk for use outside spacetime*
  }
  
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



