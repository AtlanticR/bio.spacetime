
  spacetime_db = function( DS, p, B=NULL, yr=NULL, ret="mean"  ) {
    #// usage: low level function to convert data into file-based data obects to permit parallel
    #// data access and manipulation and deletes/updates
    #// B is the xyz or xytz data or the function to get the data to work upon

    if ( DS=="spatial.covariance") {
      stats = NULL
      p = spacetime_db( p=p, DS="filenames" )
      if (file.exists( p$fn$stats) ) load( p$fn.stats )
      return(stats)
    }

    if ( DS =="predictions" ) {
      preds = NULL
      p = spacetime_db( p=p, DS="filenames" )
      if (file.exists( p$fn$P ) ) load( p$fn.P )
      return( preds )
    }
    
    if ( DS =="statistics" ) {
      stats = NULL
      p = spacetime_db( p=p, DS="filenames" )  
      if (file.exists( p$fn$S) ) load( p$fn.S )
      return( stats )
    }

    # --------------------------

    if (DS %in% "filenames" ) {
      # input data stored as a bigmemory file to permit operations with min memory usage
      # split into separate components to reduce filelocking conflicts

      p$tmp.datadir = file.path( p$project.root, "tmp" )
      # message( paste( "Temporary files are being created at:", p$tmp.datadir ) )
      if( !file.exists(p$tmp.datadir)) dir.create( p$tmp.datadir, recursive=TRUE, showWarnings=FALSE )

      p$savedir = file.path(p$project.root, "spacetime", p$spatial.domain )
      
      # message( paste( "Final outputs will be palced at:", p$savedir ) )
      if( !file.exists(p$savedir)) dir.create( p$savedir, recursive=TRUE, showWarnings=FALSE )
        
      # storage locations for finalized data 
      p$fn = list()
      p$fn$P = file.path( p$savedir, paste( "spacetime", "predictions", "rdata", sep=".") )
      p$fn$S = file.path( p$savedir, paste( "spacetime", "statistics",  "rdata", sep=".") )
      p$fn$stats =  file.path( p$project.root, "spacetime", paste( "spatial", "covariance", "rdata", sep=".") )
      p$fn$covmodel =  file.path( p$project.root, "spacetime", paste( "spatial", "covariate.model", "rdata", sep=".") )
      
      p$diskcache =list()
      p$diskcache$Y0 =    file.path( p$tmp.datadir, "input.Y0.diskcache" ) # raw data
      p$diskcache$Y =     file.path( p$tmp.datadir, "input.Y.diskcache" ) # residuals of covar model or raw data if none
      p$diskcache$Ycov =  file.path( p$tmp.datadir, "input.Ycov.diskcache"  )
      p$diskcache$Yloc =  file.path( p$tmp.datadir, "input.Yloc.diskcache" )
      p$diskcache$Ytime = file.path( p$tmp.datadir, "input.Ytime.diskcache" )
      p$diskcache$Yi =    file.path( p$tmp.datadir, "input.Yi.diskcache" ) # index of useable data
      p$diskcache$P0 =    file.path( p$tmp.datadir, "predictions0.diskcache" ) # offsets from covar model
      p$diskcache$P =     file.path( p$tmp.datadir, "predictions.diskcache" )
      p$diskcache$Psd =   file.path( p$tmp.datadir, "predictions_sd.diskcache" )
      p$diskcache$Pn =    file.path( p$tmp.datadir, "predictions_n.diskcache" )
      p$diskcache$Pcov =  file.path( p$tmp.datadir, "predictions_cov.diskcache" )
      p$diskcache$Ploc =  file.path( p$tmp.datadir, "predictions_loc.diskcache" )
      p$diskcache$Ptime = file.path( p$tmp.datadir, "predictions_time.diskcache" )
      p$diskcache$S =     file.path( p$tmp.datadir, "statistics.diskcache" )
      p$diskcache$Sloc =  file.path( p$tmp.datadir, "statistics_loc.diskcache" )
      p$diskcache$Stime = file.path( p$tmp.datadir, "statistics_time.diskcache" )
      p$diskcache$Mat2Ploc = file.path( p$tmp.datadir, "Mat2Ploc.diskcache" )
      p$diskcache$spatial_weights = file.path( p$tmp.datadir, "spatial_weights.diskcache" )

      # if ( p$storage.backend=="bigmemory.filebacked" )  {
      #   # file-backed pointers are hard coded with bigmemory
      #   p$ptr =list()
      #   p$ptr$Y0 =    file.path( p$tmp.datadir, "input.Y0.pointer" ) # raw data
      #   p$ptr$Y =     file.path( p$tmp.datadir, "input.Y.pointer" ) # residuals of covar model or raw data if none
      #   p$ptr$Ycov =  file.path( p$tmp.datadir, "input.Ycov.pointer"  )
      #   p$ptr$Yloc =  file.path( p$tmp.datadir, "input.Yloc.pointer" )
      #   p$ptr$Ytime = file.path( p$tmp.datadir, "input.Ytime.pointer" )
      #   p$ptr$Yi =    file.path( p$tmp.datadir, "input.Yi.pointer" )
      #   p$ptr$P0 =    file.path( p$tmp.datadir, "predictions0.pointer" ) # offsets from covar model
      #   p$ptr$P =     file.path( p$tmp.datadir, "predictions.pointer" )
      #   p$ptr$Psd =   file.path( p$tmp.datadir, "predictions_sd.pointer" )
      #   p$ptr$Pn =    file.path( p$tmp.datadir, "predictions_n.pointer" )
      #   p$ptr$Pcov =  file.path( p$tmp.datadir, "predictions_cov.pointer" )
      #   p$ptr$Ploc =  file.path( p$tmp.datadir, "predictions_loc.pointer" )
      #   p$ptr$Ptime = file.path( p$tmp.datadir, "predictions_time.pointer" )
      #   p$ptr$S =     file.path( p$tmp.datadir, "statistics.pointer" )
      #   p$ptr$Sloc =  file.path( p$tmp.datadir, "statistics_loc.pointer" )
      #   p$ptr$Stime = file.path( p$tmp.datadir, "statistics_time.pointer" )
      #   p$ptr$Mat2Ploc = file.path( p$tmp.datadir, "Mat2Ploc.pointer" )
      #   p$ptr$spatial_weights = file.path( p$tmp.datadir, "spatial_weights.pointer" )
      # }

      return(p)
    }

    # --------------------------

    if (DS=="save.parameters")  {
      p = spacetime_db( p=p, DS="filenames" )  
      save(p, file=file.path( p$tmp.datadir, "p.rdata") )
      message( "Saved parameters:")
      message( file.path( p$tmp.datadir, "p.rdata") )
    }
    if (DS=="load.parameters")  {
      p = spacetime_db( p=p, DS="filenames" )  
      load(file.path( p$tmp.datadir, "p.rdata") )
      RLibrary( p$libs )
      return(p)
    }

    # --------------------------
    if (DS %in% "cleanup" ) {
      p = spacetime_db( p=p, DS="filenames" )
      for (fn in p$ptr ) if (file.exists(fn)) file.remove(fn)
      for (fn in p$diskcache ) if (file.exists(fn)) file.remove(fn)
      return( "done" )
    }

    # ------------------
    if (DS == "data.initialize" ) {
      p = spacetime_db( p=p, DS="filenames" ) 
      if ( "data.frame" %in% class(B) ) {
        # B = B # nothing to do
      } else if ( "SpatialGridDataFrame" %in% class(B) ) {
        B = slot(B, "data") 
      }

      # dependent variable
      Y_ = B[, p$variables$Y ]
        if (p$storage.backend=="ff") {
          p$ptr$Y0 = ff( Y_, dim=dim(Y_), file=p$diskcache$Y0, overwrite=TRUE )
        }
        if (p$storage.backend=="bigmemory.filebacked") {
          Y0 = bigmemory::as.big.matrix( Y_, type="double", 
            backingfile=p$diskcache$Y0, descriptorfile=basename(p$ptr$Y0), backingpath=p$tmp.datadir )
          p$ptr$Y0 = bigmemory::describe(Y0) # not sure if necessary .. should check .. 
        }
        if (p$storage.backend=="bigmemory.ram") {
          Y0 = bigmemory::as.big.matrix( Y_, type="double" )
          p$ptr$Y0 = bigmemory::describe(Y0) 
        }

      if (exists( "spacetime_covariate_spacetime_engine_modelformula", p)) {
        Yresid_ = residuals( spacetime_db( p=p, DS="model.covariates") )
      } else {
        Yresid_ = Y_  # if no residual model simply use the raw data
      }
        if (p$storage.backend=="ff") {
          p$ptr$Y = ff( Yresid_, dim=dim(Yresid_), file=p$diskcache$Y, overwrite=TRUE )
        }
        if (p$storage.backend=="bigmemory.filebacked") {
          Y = bigmemory::as.big.matrix( Yresid_, type="double", 
            backingfile=p$diskcache$Y, descriptorfile=basename(p$ptr$Y), backingpath=p$tmp.datadir )
          p$ptr$Y = bigmemory::describe(Y) # not sure if necessary .. should check .. 
        }
        if (p$storage.backend=="bigmemory.ram") {
          Y = bigmemory::as.big.matrix( Yresid_, type="double" )
          p$ptr$Y = bigmemory::describe(Y) 
        }


     # data coordinates
      Yloc_ = as.matrix( B[, p$variables$LOCS ])
        if (p$storage.backend=="ff") {
          p$ptr$Yloc = ff( Yloc_, dim=dim(Yloc_), file=p$diskcache$Yloc, overwrite=TRUE )
        }
        if (p$storage.backend=="bigmemory.filebacked") {
          Yloc = bigmemory::as.big.matrix( Yloc_, type="double", 
            backingfile=p$diskcache$Yloc, descriptorfile=basename(p$ptr$Yloc), backingpath=p$tmp.datadir )
          p$ptr$Yloc = bigmemory::describe(Yloc) # not sure if necessary .. should check .. 
        }
        if (p$storage.backend=="bigmemory.ram") {
          Yloc = bigmemory::as.big.matrix( Yloc_, type="double" )
          p$ptr$Yloc = bigmemory::describe(Yloc) 
        }

      # independent variables/ covariate
      if (exists("COV", p$variables)) {
        Ycov_ = as.matrix( B[ , p$variables$COV ] )
          if (p$storage.backend=="ff") {
            p$ptr$Ycov = ff( Ycov_, dim=dim(Ycov_), file=p$diskcache$Ycov, overwrite=TRUE )
          }
          if (p$storage.backend=="bigmemory.filebacked") {
            Ycov = bigmemory::as.big.matrix( Ycov_, type="double", 
              backingfile=p$diskcache$Ycov, descriptorfile=basename(p$ptr$Ycov), backingpath=p$tmp.datadir )
            p$ptr$Ycov = bigmemory::describe(Ycov) # not sure if necessary .. should check .. 
          }
          if (p$storage.backend=="bigmemory.ram") {
            Ycov = bigmemory::as.big.matrix( Ycov_, type="double" )
            p$ptr$Ycov = bigmemory::describe(Ycov) 
          }

      }

      # data times
      if ( exists("TIME", p$variables) ) {
        Ytime_ = as.matrix( B[, p$variables$TIME ] )
          if (p$storage.backend=="ff") {
            p$ptr$Ytime = ff( Ytime_, dim=dim(Ytime_), file=p$diskcache$Ytime, overwrite=TRUE )
          }
          if (p$storage.backend=="bigmemory.filebacked") {
            Ytime = bigmemory::as.big.matrix( Ytime_, type="double", 
              backingfile=p$diskcache$Ytime, descriptorfile=basename(p$ptr$Ytime), backingpath=p$tmp.datadir )
            p$ptr$Ytime = bigmemory::describe(Ytime) # not sure if necessary .. should check .. 
          }
          if (p$storage.backend=="bigmemory.ram") {
            Ytime = bigmemory::as.big.matrix( Ytime_, type="double" )
            p$ptr$Ytime = bigmemory::describe(Ytime) 
          }
      }

      return( p ) #return pointers to data
    }

    #---------------------
    if (DS == "predictions.initialize" ) {

      p = spacetime_db( p=p, DS="filenames" ) 
      
      if (exists("COV", p$variables)) {
        if (is.vector(B$COV) ) {
          Pcov_ = as.matrix( B$COV ) 
        } else {
          Pcov_ = as.matrix( B$COV[,p$variables$COV ] ) 
        }
          if (p$storage.backend=="ff") {
            p$ptr$Pcov = ff( Pcov_, dim=dim(Pcov_), file=p$diskcache$Pcov, overwrite=TRUE )
          }
          if (p$storage.backend=="bigmemory.filebacked") {
            Pcov = bigmemory::as.big.matrix( Pcov_, type="double", 
              backingfile=p$diskcache$Pcov, descriptorfile=basename(p$ptr$Pcov), backingpath=p$tmp.datadir )
            p$ptr$Pcov = bigmemory::describe(Pcov) # not sure if necessary .. should check .. 
          }
          if (p$storage.backend=="bigmemory.ram") {
            Pcov = bigmemory::as.big.matrix( Pcov_, type="double" )
            p$ptr$Pcov = bigmemory::describe(Pcov) 
          }
      }

      # prediction times 
      if (exists("TIME", p$variables)) {
        Ptime_ = as.matrix( B$TIME )
          if (p$storage.backend=="ff") {
            p$ptr$Ptime = ff( Ptime_, dim=dim(Ptime_), file=p$diskcache$Ptime, overwrite=TRUE )
          }
          if (p$storage.backend=="bigmemory.filebacked") {
            Ptime = bigmemory::as.big.matrix( Ptime_, type="double", 
              backingfile=p$diskcache$Ptime, descriptorfile=basename(p$ptr$Ptime), backingpath=p$tmp.datadir )
            p$ptr$Ptime = bigmemory::describe(Ptime) # not sure if necessary .. should check .. 
          }
          if (p$storage.backend=="bigmemory.ram") {
            Ptime = bigmemory::as.big.matrix( Ptime_, type="double" )
            p$ptr$Ptime = bigmemory::describe(Ptime) 
          }
      }
      
      # predictions and associated stats
      P_ = matrix( NaN, nrow=nrow(B$LOCS), ncol=p$nt )
        # predictions proper
        if (p$storage.backend=="ff") {
          p$ptr$P = ff( P_, dim=dim(P_), file=p$diskcache$P, overwrite=TRUE )
        }
        if (p$storage.backend=="bigmemory.filebacked") {
          P = bigmemory::as.big.matrix( P_, type="double", 
            backingfile=p$diskcache$P, descriptorfile=basename(p$ptr$P), backingpath=p$tmp.datadir )
          p$ptr$P = bigmemory::describe(P) # not sure if necessary .. should check .. 
        }
        if (p$storage.backend=="bigmemory.ram") {
          P = bigmemory::as.big.matrix( P_, type="double" )
          p$ptr$P = bigmemory::describe(P) 
        }
        
        # count of prediction estimates
        if (p$storage.backend=="ff") {
          p$ptr$Pn = ff( P_, dim=dim(P_), file=p$diskcache$Pn, overwrite=TRUE )
        }
        if (p$storage.backend=="bigmemory.filebacked") {
          Pn = bigmemory::as.big.matrix( P_, type="double", 
            backingfile=p$diskcache$Pn, descriptorfile=basename(p$ptr$Pn), backingpath=p$tmp.datadir )
          p$ptr$Pn = bigmemory::describe(Pn) # not sure if necessary .. should check .. 
        }
        if (p$storage.backend=="bigmemory.ram") {
          Pn = bigmemory::as.big.matrix( P_, type="double" )
          p$ptr$Pn = bigmemory::describe(Pn) 
        }

        # sd of prediction estimates
        if (p$storage.backend=="ff") {
          p$ptr$Psd = ff( P_, dim=dim(P_), file=p$diskcache$Psd, overwrite=TRUE )
        }
        if (p$storage.backend=="bigmemory.filebacked") {
          Psd = bigmemory::as.big.matrix( P_, type="double", 
            backingfile=p$diskcache$Psd, descriptorfile=basename(p$ptr$Psd), backingpath=p$tmp.datadir )
          p$ptr$Psd = bigmemory::describe(Psd) # not sure if necessary .. should check .. 
        }
        if (p$storage.backend=="bigmemory.ram") {
          Psd = bigmemory::as.big.matrix( P_, type="double" )
          p$ptr$Psd = bigmemory::describe(Psd) 
        }

      # prediction coordinates
      Ploc_ = as.matrix( B$LOCS )
          if (p$storage.backend=="ff") {
            p$ptr$Ploc = ff( Ploc_, dim=dim(Ploc_), file=p$diskcache$Ploc, overwrite=TRUE )
          }
          if (p$storage.backend=="bigmemory.filebacked") {
            Ploc = bigmemory::as.big.matrix( Ploc_, type="double", 
              backingfile=p$diskcache$Ploc, descriptorfile=basename(p$ptr$Ploc), backingpath=p$tmp.datadir )
            p$ptr$Ploc = bigmemory::describe(Ploc) # not sure if necessary .. should check .. 
          }
          if (p$storage.backend=="bigmemory.ram") {
            Ploc = bigmemory::as.big.matrix( Ploc_, type="double" )
            p$ptr$Ploc = bigmemory::describe(Ploc) 
          }

      # pre-compute a few things for spacetime_interpolate_xy_simple_multiple  
      Mat2Ploc_ = cbind( (Ploc_[,1]-p$plons[1])/p$pres + 1, (Ploc_[,2]-p$plats[1])/p$pres + 1) # row, col indices in matrix form
          if (p$storage.backend=="ff") {
            p$ptr$Mat2Ploc = ff( Mat2Ploc_, dim=dim(Mat2Ploc_), file=p$diskcache$Mat2Ploc, overwrite=TRUE )
          }
          if (p$storage.backend=="bigmemory.filebacked") {
            Mat2Ploc = bigmemory::as.big.matrix( Mat2Ploc_, type="double", 
              backingfile=p$diskcache$Mat2Ploc, descriptorfile=basename(p$ptr$Mat2Ploc), backingpath=p$tmp.datadir )
            p$ptr$Mat2Ploc = bigmemory::describe(Mat2Ploc) # not sure if necessary .. should check .. 
          }
          if (p$storage.backend=="bigmemory.ram") {
            Mat2Ploc = bigmemory::as.big.matrix( Mat2Ploc_, type="double" )
            p$ptr$Mat2Ploc = bigmemory::describe(Mat2Ploc) 
          }

      spatial_weights_ = setup.image.smooth( nrow=p$nplons, ncol=p$nplats, dx=p$pres, dy=p$res, 
        theta=p$theta, xwidth=p$nsd*p$theta, ywidth=p$nsd*p$theta )
          if (p$storage.backend=="ff") {
            p$ptr$spatial_weights = ff( spatial_weights_, dim=dim(spatial_weights_), file=p$diskcache$spatial_weights, overwrite=TRUE )
          }
          if (p$storage.backend=="bigmemory.filebacked") {
            spatial_weights = bigmemory::as.big.matrix( spatial_weights_, type="double", 
              backingfile=p$diskcache$spatial_weights, descriptorfile=basename(p$ptr$spatial_weights), backingpath=p$tmp.datadir )
            p$ptr$spatial_weights = bigmemory::describe(spatial_weights) # not sure if necessary .. should check .. 
          }
          if (p$storage.backend=="bigmemory.ram") {
            spatial_weights = bigmemory::as.big.matrix( spatial_weights_, type="double" )
            p$ptr$spatial_weights = bigmemory::describe(spatial_weights) 
          }
      

      # prediction offsets from an additive (global) model of covariates (if any, zero valued otherwise) 
      if ( exists("COV", p$variables) && exists("spacetime_covariate_spacetime_engine_modelformula", p ) ) {
        # assuming model is correct..
        covmodel = spacetime_db( p=p, DS="model.covariates") 
        covars = data.frame(B$COV)
        names(covars) = p$variables$COV
        Pcov = predict( covmodel, newdata=covars, type="response", se.fit=T ) 
        P0_   = Pcov$fit
        P0sd_ = Pcov$se.fit
        rm (covmodel, covars, Pcov); gc()
      } else {
        P0_   = matrix( 0, nrow=nrow(B$LOCS), ncol=p$nt )
        P0sd_ = matrix( 0, nrow=nrow(B$LOCS), ncol=p$nt )
      }

        if (p$storage.backend=="ff") {
          p$ptr$P0 = ff( P0_, dim=dim(P0_), file=p$diskcache$P0, overwrite=TRUE )
        }
        if (p$storage.backend=="bigmemory.filebacked") {
          P0 = bigmemory::as.big.matrix( P0_, type="double", 
            backingfile=p$diskcache$P0, descriptorfile=basename(p$ptr$P0), backingpath=p$tmp.datadir )
          p$ptr$P0 = bigmemory::describe(P0) # not sure if necessary .. should check .. 
        }
        if (p$storage.backend=="bigmemory.ram") {
          P0 = bigmemory::as.big.matrix( P0_, type="double" )
          p$ptr$P0 = bigmemory::describe(P0) 
        }
      
        if (p$storage.backend=="ff") {
          p$ptr$P0sd = ff( P0sd_, dim=dim(P0sd_), file=p$diskcache$P0sd, overwrite=TRUE )
        }
        if (p$storage.backend=="bigmemory.filebacked") {
          P0sd = bigmemory::as.big.matrix( P0sd_, type="double", 
            backingfile=p$diskcache$P0sd, descriptorfile=basename(p$ptr$P0sd), backingpath=p$tmp.datadir )
          p$ptr$P0sd = bigmemory::describe(P0sd) # not sure if necessary .. should check .. 
        }
        if (p$storage.backend=="bigmemory.ram") {
          P0sd = bigmemory::as.big.matrix( P0sd_, type="double" )
          p$ptr$P0sd = bigmemory::describe(P0sd) 
        }
      
      return( p )
    }

    # -----------------

    if (DS %in% c( "boundary.redo", "boundary" ) )  {
      p = spacetime_db( p=p, DS="filenames" )
      fn =  file.path(p$tmp.datadir, "boundary.rdata" )
      if (DS=="boundary") {
        boundary = NULL
        if( file.exists(fn)) load( fn)
        return( boundary )
      }

      # data:
      Y = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$Y), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$Y), 
        ff=p$ptr$Y )
      hasdata = 1:length(Y)
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA

      # covariates (independent vars)
      if ( exists( "COV", p$variables) ) {
        Ycov = switch( p$storage.backend, 
          bigmemory.ram=attach.big.matrix(p$ptr$Ycov), 
          bigmemory.filebacked=attach.big.matrix(p$ptr$Ycov), 
          ff=p$ptr$Ycov )
 
        if ( length( p$variables$COV ) == 1 ) {
          bad = which( !is.finite( Ycov[]) )
        } else {
          bad = which( !is.finite( rowSums(Ycov[])) )
        }
        if (length(bad)> 0 ) hasdata[bad] = NA
      }

      ii = na.omit(hasdata)
      ndata = length(ii)
      Yloc = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$Yloc), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$Yloc), 
        ff=p$ptr$Yloc )

      locs_noise =  runif( ndata*2, min=-p$pres*p$spacetime.noise, max=p$pres*p$spacetime.noise )
      if (!exists("non_convex_hull_alpha", p)) p$non_convex_hull_alpha=20
      boundary=list( polygon = non_convex_hull( Yloc[ii,]+locs_noise, alpha=p$non_convex_hull_alpha, plot=FALSE ) )
      
      # statistical output locations
      Sloc = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$Sloc), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$Sloc), 
        ff=p$ptr$Sloc )
 
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

    # -----------------
    
    if (DS %in% c( "statistics.initialize", "statistics.status", "statistics.box" ) ) {
      p = spacetime_db( p=p, DS="filenames" )  
     
      if ( DS=="statistics.initialize" ) {
        # statistics storage matrix ( aggregation window, coords ) .. no inputs required
        # statistics coordinates
        Sloc_ = as.matrix( expand.grid( p$sbbox$plons, p$sbbox$plats ))
          if (p$storage.backend=="ff") {
            p$ptr$Sloc = ff( Sloc_, dim=dim(Sloc_), filename=p$diskcache$Sloc, overwrite=TRUE )
          }
          if (p$storage.backend=="bigmemory.filebacked") {
            Sloc = bigmemory::as.big.matrix( Sloc_, type="double", 
              backingfile=p$diskcache$Sloc, descriptorfile=basename(p$ptr$Sloc), backingpath=p$tmp.datadir )
            p$ptr$Sloc = bigmemory::describe(Sloc) # not sure if necessary .. should check .. 
          }
          if (p$storage.backend=="bigmemory.ram") {
            Sloc = bigmemory::as.big.matrix( Sloc_, type="double" )
            p$ptr$Sloc = bigmemory::describe(Sloc) 
          }

        S_ = matrix( NaN, nrow=nrow(Sloc_), ncol=length( p$statsvars ) ) # NA forces into logical
          if (p$storage.backend=="ff") {
            p$ptr$S = ff( S_, dim=dim(S_), file=p$diskcache$S, overwrite=TRUE )
          }
          if (p$storage.backend=="bigmemory.filebacked") {
            S = bigmemory::as.big.matrix( S_, type="double", 
              backingfile=p$diskcache$S, descriptorfile=basename(p$ptr$S), backingpath=p$tmp.datadir )
            p$ptr$S = bigmemory::describe(S) # not sure if necessary .. should check .. 
          }
          if (p$storage.backend=="bigmemory.ram") {
            S = bigmemory::as.big.matrix( S_, type="double" )
            p$ptr$S = bigmemory::describe(S) 
          }

        return( p )
      }
      
      if (DS == "statistics.box")  {
        sbbox = list( plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$spacetime.prediction.dist.min ),
                      plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$spacetime.prediction.dist.min ) )
        return(sbbox)
      }


      if ( DS=="statistics.status" ) {
        # find locations for statistic computation and trim area based on availability of data
        # stats:
        S = switch( p$storage.backend, 
          bigmemory.ram=attach.big.matrix(p$ptr$S), 
          bigmemory.filebacked=attach.big.matrix(p$ptr$S), 
          ff=p$ptr$S )

        i = which( is.infinite( S[,1] )  )  # not yet completed (due to a failed attempt)
        j = which( is.nan( S[,1] )   )      # completed
        k = which( is.finite (S[,1])  )     # not yet done
        out = list(problematic=i, incomplete=j, completed=k, n.total=nrow(S) ,
                     n.incomplete=length(j), n.problematic=length(i), 
                     n.complete=length(k) )
        out$prop_incomp=out$n.incomplete / ( out$n.incomplete + out$n.complete)
        message( paste("Proportion incomplete:", round(out$prop_incomp,5), "\n" )) 
        return( out )
      }
    }

    # -----

    if (DS=="data.filter") {
      # last set of filters to reduce problem size
      S = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$S), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$S), 
        ff=p$ptr$S )

      bnds = try( spacetime_db( p, DS="boundary" ) )
      if (!is.null(bnds)) {
        if( !("try-error" %in% class(bnds) ) ) {
          # problematic and/or no data (e.g., land, etc.) and skipped
          to.ignore =  which( bnds$inside.polygon == 0 ) # outside boundary
          if (length(to.ignore)>0) S[to.ignore,] = Inf
      }}

      Y = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$Y), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$Y), 
        ff=p$ptr$Y )

      Yloc = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$Yloc), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$Yloc), 
        ff=p$ptr$Yloc )

      Yi = 1:length(Y) # index with useable data
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) Yi[bad] = NA

      # data locations
      bad = which( !is.finite( rowSums(Yloc[])))
      if (length(bad)> 0 ) Yi[bad] = NA

    # data locations
      if (exists("COV", p$variables)) {
        Ycov = switch( p$storage.backend, 
          bigmemory.ram=attach.big.matrix(p$ptr$Ycov), 
          bigmemory.filebacked=attach.big.matrix(p$ptr$Ycov), 
          ff=p$ptr$Ycov )

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
        Ytime = switch( p$storage.backend, 
          bigmemory.ram=attach.big.matrix(p$ptr$Ytime), 
          bigmemory.filebacked=attach.big.matrix(p$ptr$Ytime), 
          ff=p$ptr$Ytime )
        bad = which( !is.finite( Ytime[] ))
        if (length(bad)> 0 ) Yi[bad] = NA
        Yi = na.omit(Yi)
      }

          if (p$storage.backend=="ff") {
            p$ptr$Yi = ff( Yi, filename=p$diskcache$Yi, overwrite=TRUE )
          }
          if (p$storage.backend=="bigmemory.filebacked") {
            Yi = bigmemory::as.big.matrix( Yi, type="double", 
              backingfile=p$diskcache$Yi, descriptorfile=basename(p$ptr$Yi), backingpath=p$tmp.datadir )
            p$ptr$Yi = bigmemory::describe(Yi) # not sure if necessary .. should check .. 
          }
          if (p$storage.backend=="bigmemory.ram") {
            Yi = bigmemory::as.big.matrix( Yi, type="double" )
            p$ptr$Yi = bigmemory::describe(Yi) 
          }


      #---------------------
      # prediction locations and covariates
      Ploc = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$Ploc), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$Ploc), 
        ff=p$ptr$Ploc )

      rcP = data.frame( cbind( 
        Prow = (Ploc[,1]-p$plons[1])/p$pres + 1,  
        Pcol = (Ploc[,2]-p$plats[1])/p$pres + 1) )
      rcP$rc = paste( rcP$Prow, rcP$Pcol, sep="~")
      rcP$Prow = rcP$Pcol = NULL
      p$rcP = rcP

      #-----------------
      # row, col indices
      # statistical output locations
      Sloc = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$Sloc), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$Sloc), 
        ff=p$ptr$Sloc )
      
      rcS = data.frame( cbind( 
        Srow = (Sloc[,1]-p$plons[1])/p$pres + 1,  
        Scol = (Sloc[,2]-p$plats[1])/p$pres + 1))
      p$rcS = rcS
      return(p)

    }


    # -----
  
    if (DS %in% c("model.covariates", "model.covariates.redo") ) {
      
      if (DS =="model.covariates") {
        covmodel = NULL
        if (file.exists( p$fn$covmodel ))  load(p$fn$covmodel)
        return(covmodel)
      }  
      
      # as a first pass, model the time-independent factors as a user-defined model
      if (p$spacetime_covariate_modeltype=="gam") {
        covmodel = try( 
          gam( p$spacetime_covariate_spacetime_engine_modelformula, data=B, optimizer=c("outer","bfgs") ) ) 
        if ( "try-error" %in% class(covmodel) ) stop( "The covariate model was problematic" )
        print( summary( covmodel ) )
        save( covmodel, file= p$fn$covmodel, compress=TRUE )
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

    
      PP = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$P), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$P), 
        ff=p$ptr$P )

      PPsd = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$Psd), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$Psd), 
        ff=p$ptr$Psd )

      P0 = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$P0), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$P0), 
        ff=p$ptr$P0 )

      P0sd = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$P0sd), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$P0sd), 
        ff=p$ptr$P0sd )

      PP = P0 + PP 
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

    if (DS %in% c("stats_to_prediction_grid.redo", "stats_to_prediction_grid") ) {

      if (DS=="stats_to_prediction_grid") {
        stats = NULL
        if (file.exists(p$fn.S)) {}
      }

      S = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$S), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$S), 
        ff=p$ptr$S )

      Sloc = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$Sloc), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$Sloc), 
        ff=p$ptr$Sloc )
    
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
          theta=p$spacetime.prediction.dist.min, xwidth=p$spacetime.prediction.dist.min*10, ywidth=p$spacetime.prediction.dist.min*10 )
        if (!is.null(res)) stats[i,] = res
      }
    
      # subset to match to Ploc
      locsout_rc = paste( locsout$plon, locsout$plat, sep="~" )
      Ploc = switch( p$storage.backend, 
        bigmemory.ram=attach.big.matrix(p$ptr$Ploc), 
        bigmemory.filebacked=attach.big.matrix(p$ptr$Ploc), 
        ff=p$ptr$Ploc )

      bad = which( !is.finite(pa$i))
      if (length(bad) > 0 ) pa = pa[-bad,]
      if (nrow(pa)< 5) next()
      pa$plon = Ploc[ pa$i, 1]
      pa$plat = Ploc[ pa$i, 2]

      save( stats,  file=p$fn.S, compress=TRUE )

    }

  }


