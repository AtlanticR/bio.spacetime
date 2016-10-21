
  spacetime_db = function( DS, p, B=NULL, yr=NULL, ret="mean"  ) {
    #// usage: low level function to convert data into file-based data obects to permit parallel
    #// data access and manipulation and deletes/updates
    #// B is the xyz or xytz data or the function to get the data to work upon

    # --------------------------

    if (DS %in% "filenames" ) {
      # input data stored as a bigmemory file to permit operations with min memory usage
      # split into separate components to reduce filelocking conflicts

      p$stloc = file.path( p$project.root, "tmp" )
      # message( paste( "Temporary files are being created at:", p$stloc ) )
      if( !file.exists(p$stloc)) dir.create( p$stloc, recursive=TRUE, showWarnings=FALSE )

      p$savedir = file.path(p$project.root, "spacetime", p$spatial.domain )
      
      # message( paste( "Final outputs will be palced at:", p$savedir ) )
      if( !file.exists(p$savedir)) dir.create( p$savedir, recursive=TRUE, showWarnings=FALSE )
        
      # storage locations for finalized data 
      p$fn = list()
      p$fn$P = file.path( p$savedir, paste( "spacetime", "predictions", "rdata", sep=".") )
      p$fn$S = file.path( p$savedir, paste( "spacetime", "statistics",  "rdata", sep=".") )
      p$fn$stats =  file.path( p$project.root, "spacetime", paste( "spatial", "covariance", "rdata", sep=".") )
     
      p$cache =list()
      p$cache$Y0 =    file.path( p$stloc, "input.Y0.cache" ) # raw data
      p$cache$Y =     file.path( p$stloc, "input.Y.cache" ) # residuals of covar model or raw data if none
      p$cache$Ycov =  file.path( p$stloc, "input.Ycov.cache"  )
      p$cache$Yloc =  file.path( p$stloc, "input.Yloc.cache" )
      p$cache$Ytime = file.path( p$stloc, "input.Ytime.cache" )
      p$cache$Yi =    file.path( p$stloc, "input.Yi.cache" ) # index of useable data
      p$cache$P0 =    file.path( p$stloc, "predictions0.cache" ) # offsets from covar model
      p$cache$P0sd =  file.path( p$stloc, "predictions0sd.cache" ) # offsets from covar model
      p$cache$P =     file.path( p$stloc, "predictions.cache" )
      p$cache$Psd =   file.path( p$stloc, "predictions_sd.cache" )
      p$cache$Pn =    file.path( p$stloc, "predictions_n.cache" )
      p$cache$Pcov =  file.path( p$stloc, "predictions_cov.cache" )
      p$cache$Ploc =  file.path( p$stloc, "predictions_loc.cache" )
      p$cache$Ptime = file.path( p$stloc, "predictions_time.cache" )
      p$cache$S =     file.path( p$stloc, "statistics.cache" )
      p$cache$Sloc =  file.path( p$stloc, "statistics_loc.cache" )
      p$cache$Stime = file.path( p$stloc, "statistics_time.cache" )
      p$cache$Mat2Ploc = file.path( p$stloc, "Mat2Ploc.cache" )
      p$cache$Ylogit = file.path( p$stloc, "Ylogit.cache" )

      p$bm = p$cache
      for (i in 1:length(p$bm)) p$bm[[i]] = gsub(".cache$", ".bigmemory", p$bm[[i]] )

      return(p)
    }

    # --------------------------

    if (DS=="save.parameters")  {
      save(p, file=file.path( p$stloc, "p.rdata") )
      message( "Saved parameters:")
      message( file.path( p$stloc, "p.rdata") )
    }
    if (DS=="load.parameters")  {
      load(file.path( p$stloc, "p.rdata") )
      RLibrary( p$libs )
      return(p)
    }

    # --------------------------
    if (DS %in% "cleanup" ) {
      for (fn in p$cache ) if (file.exists(fn)) file.remove(fn)
      for (fn in p$bm ) if (file.exists(fn)) file.remove(fn)
      return( "done" )
    }

    # -----------------
    
    if (DS %in% c( "statistics.initialize", "statistics.status", "statistics.box" ) ) {
     
      if ( DS=="statistics.initialize" ) {
        # statistics storage matrix ( aggregation window, coords ) .. no inputs required
        # statistics coordinates
        Sloc = as.matrix( expand.grid( p$sbbox$plons, p$sbbox$plats ))
          if (p$storage.backend == "bigmemory.ram" ) {
            p$ptr$Sloc  = bigmemory::describe( bigmemory::as.big.matrix( Sloc, type="double" ) )
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
            p$ptr$S  = bigmemory::describe( bigmemory::as.big.matrix( S, type="double" ) )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$S  = p$cache$S
            bigmemory::as.big.matrix( S, type="double", backingfile=basename(p$bm$S), descriptorfile=basename(p$cache$S), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$S = ff( S, dim=dim(S), file=p$cache$S, overwrite=TRUE )
          }

        return( p )
      }
      
      if (DS == "statistics.box")  {
        sbbox = list( plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$spacetime_prediction_dist_min ),
                      plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$spacetime_prediction_dist_min ) )
        return(sbbox)
      }


      if ( DS=="statistics.status" ) {
        # find locations for statistic computation and trim area based on availability of data
        # stats:
        S = spacetime_attach( p$storage.backend, p$ptr$S )
        i = which( is.infinite( S[,1] )  )  # not yet completed (due to a failed attempt)
        j = which( is.nan( S[,1] )   )      # incomplete
        k = which( is.finite (S[,1])  )     # completed
        out = list(problematic=i, incomplete=j, completed=k, n.total=nrow(S) ,
                     n.incomplete=length(j), n.problematic=length(i), 
                     n.complete=length(k) )
        out$prop_incomp=out$n.incomplete / ( out$n.incomplete + out$n.complete)
        message( paste("Proportion incomplete:", round(out$prop_incomp,5), "\n" )) 
        return( out )
      }
    }

    # ------------------
    if (DS == "data.initialize" ) {
      if ( "data.frame" %in% class(B) ) {
        # B = B # nothing to do
      } else if ( "SpatialGridDataFrame" %in% class(B) ) {
        B = slot(B, "data") 
      }

      # dependent variable
      Y = as.matrix(B[, p$variables$Y ])
        if (p$storage.backend == "bigmemory.ram" ) {
          p$ptr$Y  = bigmemory::describe( bigmemory::as.big.matrix( Y, type="double" ) )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Y  = p$cache$Y
          bigmemory::as.big.matrix( Y, type="double", backingfile=basename(p$bm$Y), descriptorfile=basename(p$cache$Y), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Y = ff( Y, dim=dim(Y), file=p$cache$Y, overwrite=TRUE )
        }
        

      if (p$spacetime_method=="habitat") {
        if (p$storage.backend == "bigmemory.ram" ) {
          p$ptr$Ylogit  = bigmemory::describe( bigmemory::as.big.matrix( Ylogit, type="double" ) )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Ylogit  = p$cache$Ylogit
          bigmemory::as.big.matrix( Ylogit, type="double", backingfile=basename(p$bm$Ylogit), descriptorfile=basename(p$cache$Ylogit), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Ylogit = ff( Ylogit, dim=dim(Ylogit), file=p$cache$Ylogit, overwrite=TRUE )
        }
      }

     # data coordinates
      Yloc = as.matrix( B[, p$variables$LOCS ])
        if (p$storage.backend == "bigmemory.ram" ) {
          p$ptr$Yloc  = bigmemory::describe(  bigmemory::as.big.matrix( Yloc, type="double" ) )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Yloc  = p$cache$Yloc
          bigmemory::as.big.matrix( Yloc, type="double", backingfile=basename(p$bm$Yloc), descriptorfile=basename(p$cache$Yloc), backingpath=p$stloc )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Yloc = ff( Yloc, dim=dim(Yloc), file=p$cache$Yloc, overwrite=TRUE )
        }
        
      # independent variables/ covariate
      if (exists("COV", p$variables)) {
        Ycov = as.matrix( B[ , p$variables$COV ] )
          if (p$storage.backend == "bigmemory.ram" ) {
            p$ptr$Ycov  = bigmemory::describe(  bigmemory::as.big.matrix( Ycov, type="double" ) )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ycov  = p$cache$Ycov
            bigmemory::as.big.matrix( Ycov, type="double", backingfile=basename(p$bm$Ycov), descriptorfile=basename(p$cache$Ycov), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ycov = ff( Ycov, dim=dim(Ycov), file=p$cache$Ycov, overwrite=TRUE )
          }
      }

      # data times
      if ( exists("TIME", p$variables) ) {
        Ytime = as.matrix( B[, p$variables$TIME ] )
          if (p$storage.backend == "bigmemory.ram" ) {
            p$ptr$Ytime  = bigmemory::describe( bigmemory::as.big.matrix( Ytime, type="double" ) )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ytime  = p$cache$Ytime
            bigmemory::as.big.matrix( Ytime, type="double", backingfile=basename(p$bm$Ytime), descriptorfile=basename(p$cache$Ytime), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ytime = ff( Ytime, dim=dim(Ytime), file=p$cache$Ytime, overwrite=TRUE )
          }
      }

      return( p ) #return pointers to data
    }

    #---------------------
    if (DS == "predictions.initialize" ) {

      if (exists("COV", p$variables)) {
        if (is.vector(B$COV) ) {
          Pcov = as.matrix( B$COV ) 
        } else {
          Pcov = as.matrix( B$COV[,p$variables$COV ] ) 
        }
        attr( Pcov, "dimnames" ) = NULL
          if (p$storage.backend == "bigmemory.ram" ) {
            p$ptr$Pcov  = bigmemory::describe( bigmemory::as.big.matrix( Pcov, type="double" ) )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Pcov  = p$cache$Pcov
            bigmemory::as.big.matrix( Pcov, type="double", backingfile=basename(p$bm$Pcov), descriptorfile=basename(p$cache$Pcov), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Pcov = ff( Pcov, dim=dim(Pcov), file=p$cache$Pcov, overwrite=TRUE )
          }
      }
      rm(Pcov)

      # prediction times 
      if (exists("TIME", p$variables)) {
        Ptime = as.matrix( B$TIME )
        attr( Ptime, "dimnames" ) = NULL
          if (p$storage.backend == "bigmemory.ram" ) {
            p$ptr$Ptime  = bigmemory::describe( bigmemory::as.big.matrix( Ptime, type="double" ) )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ptime  = p$cache$Ptime
            bigmemory::as.big.matrix( Ptime, type="double", backingfile=basename(p$bm$Ptime), descriptorfile=basename(p$cache$Ptime), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ptime = ff( Ptime, dim=dim(Ptime), file=p$cache$Ptime, overwrite=TRUE )
          }
      }
      rm(Ptime)

      # predictions and associated stats
      P = matrix( NaN, nrow=nrow(B$LOCS), ncol=p$nt )
          if (p$storage.backend == "bigmemory.ram" ) {
            p$ptr$P  = bigmemory::describe( bigmemory::as.big.matrix( P, type="double" ) )
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
            p$ptr$Pn = bigmemory::describe( bigmemory::as.big.matrix( P, type="double" ) )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Pn  = p$cache$Pn
            bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Pn), descriptorfile=basename(p$cache$Pn), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Pn = ff( P, dim=dim(P), file=p$cache$Pn, overwrite=TRUE )
          }

      # sd of prediction estimates
      # count of prediction estimates
          if (p$storage.backend == "bigmemory.ram" ) {
            p$ptr$Psd =bigmemory::describe( bigmemory::as.big.matrix( P, type="double" ) )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Psd  = p$cache$Psd
            bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Psd), descriptorfile=basename(p$cache$Psd), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Psd = ff( P, dim=dim(P), file=p$cache$Psd, overwrite=TRUE )
          }

      rm(P)

      # prediction coordinates
      Ploc = as.matrix( B$LOCS )
      attr( Ploc, "dimnames" ) = NULL
         if (p$storage.backend == "bigmemory.ram" ) {
            p$ptr$Ploc  = bigmemory::describe( bigmemory::as.big.matrix( Ploc, type="double" ) )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Ploc  = p$cache$Ploc
            bigmemory::as.big.matrix( Ploc, type="double", backingfile=basename(p$bm$Ploc), descriptorfile=basename(p$cache$Ploc), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ploc = ff( Ploc, dim=dim(Ploc), file=p$cache$Ploc, overwrite=TRUE )
          }

      # pre-compute a few things for spacetime_interpolate_xy_simple_multiple  
      Mat2Ploc = as.matrix( cbind( (Ploc[,1]-p$plons[1])/p$pres + 1, (Ploc[,2]-p$plats[1])/p$pres + 1) ) # row, col indices in matrix form
      attr( Mat2Ploc, "dimnames" ) = NULL
          if (p$storage.backend == "bigmemory.ram" ) {
            p$ptr$Mat2Ploc  = bigmemory::describe( bigmemory::as.big.matrix( Mat2Ploc, type="double" ) )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Mat2Ploc  = p$cache$Mat2Ploc
            bigmemory::as.big.matrix( Mat2Ploc, type="double", backingfile=basename(p$bm$Mat2Ploc), descriptorfile=basename(p$cache$Mat2Ploc), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Mat2Ploc = ff( Mat2Ploc, dim=dim(Mat2Ploc), file=p$cache$Mat2Ploc, overwrite=TRUE )
          }

      rm(Mat2Ploc, Ploc)

      p$spatial_weights = setup.image.smooth( nrow=p$nplons, ncol=p$nplats, dx=p$pres, dy=p$pres, 
        theta=p$theta, xwidth=p$nsd*p$theta, ywidth=p$nsd*p$theta )

      P0   = matrix( 0, nrow=nrow(B$LOCS), ncol=p$nt )
          if (p$storage.backend == "bigmemory.ram" ) {
            p$ptr$P0  = bigmemory::describe( bigmemory::as.big.matrix( P0, type="double" ) )
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
            p$ptr$P0sd  = bigmemory::describe( bigmemory::as.big.matrix( P0, type="double" ) )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$P0sd  = p$cache$P0sd
            bigmemory::as.big.matrix( P0, type="double", backingfile=basename(p$bm$P0sd), descriptorfile=basename(p$cache$P0sd), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$P0sd = ff( P0, dim=dim(P0), file=p$cache$P0sd, overwrite=TRUE )
          }

      gc()

      return( p )
    }

    # -----------------

    if (DS %in% c( "boundary.redo", "boundary" ) )  {

      fn =  file.path(p$stloc, "boundary.rdata" )
      if (DS=="boundary") {
        boundary = NULL
        if( file.exists(fn)) load( fn)
        return( boundary )
      }

      # data:
      Y = spacetime_attach(  p$storage.backend, p$ptr$Y )
      hasdata = 1:length(Y)
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA

      # covariates (independent vars)
      if ( exists( "COV", p$variables) ) {
        Ycov = spacetime_attach(  p$storage.backend, p$ptr$Ycov )
        if ( length( p$variables$COV ) == 1 ) {
          bad = which( !is.finite( Ycov[]) )
        } else {
          bad = which( !is.finite( rowSums(Ycov[])) )
        }
        if (length(bad)> 0 ) hasdata[bad] = NA
      }

      ii = na.omit(hasdata)
      ndata = length(ii)
      Yloc = spacetime_attach(  p$storage.backend, p$ptr$Yloc )
      locs_noise =  runif( ndata*2, min=-p$pres*p$spacetime.noise, max=p$pres*p$spacetime.noise )
      if (!exists("non_convex_hull_alpha", p)) p$non_convex_hull_alpha=20
      boundary=list( polygon = non_convex_hull( Yloc[ii,]+locs_noise, alpha=p$non_convex_hull_alpha, plot=FALSE ) )
      
      # statistical output locations
      Sloc = spacetime_attach(  p$storage.backend, p$ptr$Sloc )
    
      boundary$inside.polygon = point.in.polygon( Sloc[,1], Sloc[,2],
          boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )
      
      save( boundary, file=fn, compress=TRUE )
      plot( Yloc[ii,], pch=".", col="grey" ) # data locations
      points( Sloc[which(boundary$inside.polygon==1),], pch=".", col="orange" )
      lines( boundary$polygon[] , col="green", pch=2 )
      message( "Check the map of data and boundaries. ")
      message( "If not suitable, set another value for p$non_convex_hull_alpha value (radius; distance) ")
      message( "and re-run spacetime() with the flag: overwrite=FALSE " )
      return( fn )
    }


    # -----

    if (DS=="data.filter") {
      # last set of filters to reduce problem size
      S = spacetime_attach( p$storage.backend, p$ptr$S )
      bnds = try( spacetime_db( p, DS="boundary" ) )
      if (!is.null(bnds)) {
        if( !("try-error" %in% class(bnds) ) ) {
          # problematic and/or no data (e.g., land, etc.) and skipped
          to.ignore =  which( bnds$inside.polygon == 0 ) # outside boundary
          if (length(to.ignore)>0) S[to.ignore,] = Inf
      }}

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

      Yi = as.matrix(Yi)
          if (p$storage.backend == "bigmemory.ram" ) {
            p$ptr$Yi  = bigmemory::describe( bigmemory::as.big.matrix( Yi, type="double" ) )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Yi  = p$cache$Yi
            bigmemory::as.big.matrix( Yi, type="double", backingfile=basename(p$bm$Yi), descriptorfile=basename(p$cache$Yi), backingpath=p$stloc )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Yi = ff( Yi, dim=dim(Yi), file=p$cache$Yi, overwrite=TRUE )
          }


      #---------------------
      # prediction locations and covariates
      Ploc = spacetime_attach( p$storage.backend, p$ptr$Ploc )
      rcP = data.frame( cbind( 
        Prow = (Ploc[,1]-p$plons[1])/p$pres + 1,  
        Pcol = (Ploc[,2]-p$plats[1])/p$pres + 1) )
      rcP$rc = paste( rcP$Prow, rcP$Pcol, sep="~")
      rcP$Prow = rcP$Pcol = NULL
      p$rcP = rcP

      #-----------------
      # row, col indices
      # statistical output locations
      Sloc = spacetime_attach( p$storage.backend, p$ptr$Sloc )
      rcS = data.frame( cbind( 
        Srow = (Sloc[,1]-p$plons[1])/p$pres + 1,  
        Scol = (Sloc[,2]-p$plats[1])/p$pres + 1))
      p$rcS = rcS
      return(p)

    }

    # -----
  
    if (DS %in% c("model.covariates", "model.covariates.redo") ) {
      
      fn.covmodel =  file.path( p$project.root, "spacetime", paste( "spatial", "covariate.model", p$spacetime_method, p$spacetime_family, "rdata", sep=".") )

      if (DS =="model.covariates") {
        covmodel = NULL
        if (file.exists( fn.covmodel ))  load(fn.covmodel)
        return(covmodel)
      }  
   
      # as a first pass, model the time-independent factors as a user-defined model
      if (p$spacetime_covariate_modeltype=="gam") {
        covmodel = try( 
          gam( p$spacetime_covariate_modelformula, data=B, optimizer=c("outer","bfgs"), family=p$spacetime_family ) ) 
        if ( "try-error" %in% class(covmodel) ) stop( "The covariate model was problematic" )
        print( summary( covmodel ) )
        save( covmodel, file= fn.covmodel, compress=TRUE )
      }
     
      return ( p )
    }

    # -----
    
    if (DS %in% c("spacetime.prediction.redo", "spacetime.prediction") )  {

      if (DS=="spacetime.prediction")  {
        if (! exists("TIME", p)) yr = "0000"
        if (!is.null(yr) ) {
          if (ret=="mean") {
            fn = file.path( p$savedir, paste("spacetime.prediction.mean", yr, "rdata", sep="." ) )
            if (file.exists(fn) ) load(fn)
            return (P)
          }
          if (ret=="sd") {
            fn = file.path( p$savedir, paste("spacetime.prediction.sd", yr, "rdata", sep="." ) )
            if (file.exists(fn) ) load( fn )
            return( Psd)
          }
        }
      }

      PP = spacetime_attach( p$storage.backend, p$ptr$P )
      PPsd = spacetime_attach( p$storage.backend, p$ptr$Psd )
      P0 = spacetime_attach( p$storage.backend, p$ptr$P0 )
      P0sd = spacetime_attach( p$storage.backend, p$ptr$P0sd )
    
      PP = PP + P0 
      PPsd = sqrt( P0sd^2 + PPsd^2) # simpleadditive independent errors assumed

      if ( exists("TIME", p)) {
        for ( r in 1:length(p$tyears) ) {
          y = p$tyears[r]
          fn1 = file.path( p$savedir, paste("spacetime.prediction.mean",  y, "rdata", sep="." ) )
          fn2 = file.path( p$savedir, paste("spacetime.prediction.sd",  y, "rdata", sep="." ) )
          if (exists("nw", p)) {
            col.ranges = (r-1) * p$nw + (1:p$nw) 
            P = PP  [,col.ranges]
            V = PPsd[,col.ranges] # simpleadditive independent errors assumed
          } else {
            P = PP[,r]
            V = PPsd[,r]
          }
          save( P, file=fn1, compress=T )
          save( V, file=fn2, compress=T )
          print ( paste("Year:", y)  )
        } 
      } else {
          y = "0000"
          fn1 = file.path( p$savedir, paste("spacetime.prediction.mean",  y, "rdata", sep="." ) )
          fn2 = file.path( p$savedir, paste("spacetime.prediction.sd",  y, "rdata", sep="." ) )
          save( P, file=fn1, compress=T )
          save( V, file=fn2, compress=T )
      }
    }
  
    # ----------------

    if (DS %in% c("stats.to.prediction.grid.redo", "stats.to.prediction.grid") ) {

      if (DS=="stats.to.prediction.grid") {
        stats = NULL
        if (file.exists(p$fn.S)) {}
      }

      S = spacetime_attach( p$storage.backend, p$ptr$S )
      Sloc = spacetime_attach( p$storage.backend, p$ptr$Sloc )
    
      ss = as.data.frame( cbind( Sloc[], S[] ) )
      names(ss) = c( p$variables$LOCS, p$statsvars )
      locsout = expand.grid( p$plons, p$plats ) # final output grid
      attr( locsout , "out.attrs") = NULL
      names( locsout ) = p$variables$LOCS
      stats = matrix( NaN, ncol=nstatvars, nrow=nrow( locsout) )  # output data
      colnames(stats)=p$statsvars
    
      for ( i in 1:nstatvars ) {
        data = list( x=p$sbbox$plons, y=p$sbbox$plats, z=S[,i] )
        res = spacetime_interpolate_xy_simple( interp.method="kernel.density", 
          data=ss, locsout=locsout, nr=length(p$plons), nc=length( p$plats),  
          theta=p$spacetime_prediction_dist_min, xwidth=p$spacetime_prediction_dist_min*10, ywidth=p$spacetime_prediction_dist_min*10 )
        if (!is.null(res)) stats[i,] = res
      }
    
      # subset to match to Ploc
      locsout_rc = paste( locsout$plon, locsout$plat, sep="~" )
      Ploc = spacetime_attach( p$storage.backend, p$ptr$Ploc )

      bad = which( !is.finite(pa$i))
      if (length(bad) > 0 ) pa = pa[-bad,]
      if (nrow(pa)< 5) next()
      pa$plon = Ploc[ pa$i, 1]
      pa$plat = Ploc[ pa$i, 2]

      save( stats,  file=p$fn.S, compress=TRUE )

    }

  }


