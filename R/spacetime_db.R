
  spacetime_db = function( DS, p, B=NULL, grp=NULL, yr=NULL ) {
    #// usage: low level function to convert data into file-based data obects to permit parallel
    #// data access and manipulation and deletes/updates
    #// B is the xyz or xytz data or the function to get the data to work upon

    if ( DS=="spatial.covariance") {
      stats = NULL
      p = spacetime_db( p=p, DS="filenames" )
      if (file.exists( p$fn.stats) ) load( p$fn.stats )
      return(stats)
    }
    if ( DS =="predictions" ) {
      preds = NULL
      p = spacetime_db( p=p, DS="filenames" )
      if (file.exists( p$fn.P ) ) load( p$fn.P )
      return( preds )
    }
    if ( DS =="statistics" ) {
      stats = NULL
      p = spacetime_db( p=p, DS="filenames" )  
      if (file.exists( p$fn.S) ) load( p$fn.S )
      return( stats )
    }

    # --------------------------

    if (DS %in% "filenames" ) {
      # input data stored as a bigmemory file to permit operations with min memory usage
      # split into separate components to reduce filelocking conflicts

      # storage locations for finalized data 
      p$fn = list()
      p$fn$P = file.path( p$savedir, paste( "spacetime", "predictions", "rdata", sep=".") )
      p$fn$S = file.path( p$savedir, paste( "spacetime", "statistics",  "rdata", sep=".") )
      p$fn$results.covar =  file.path( p$project.root, "spacetime", paste( "spatial", "covariance", "rdata", sep=".") )
      
      p$ptr =list()
      p$ptr$Y =     file.path( p$tmp.datadir, "input.Y.bigmemory_description" )
      p$ptr$Ycov =  file.path( p$tmp.datadir, "input.Ycov.bigmemory_description"  )
      p$ptr$Yloc =  file.path( p$tmp.datadir, "input.Yloc.bigmemory_description" )
      p$ptr$Ytime = file.path( p$tmp.datadir, "input.Ytime.bigmemory_description" )
      p$ptr$P =     file.path( p$tmp.datadir, "predictions.bigmemory_description" )
      p$ptr$Psd =   file.path( p$tmp.datadir, "predictions_sd.bigmemory_description" )
      p$ptr$Pn =    file.path( p$tmp.datadir, "predictions_n.bigmemory_description" )
      p$ptr$Pcov =  file.path( p$tmp.datadir, "predictions_cov.bigmemory_description" )
      p$ptr$Ploc =  file.path( p$tmp.datadir, "predictions_loc.bigmemory_description" )
      p$ptr$Ptime = file.path( p$tmp.datadir, "predictions_time.bigmemory_description" )
      p$ptr$S =     file.path( p$tmp.datadir, "statistics.bigmemory_description" )
      p$ptr$Sloc =  file.path( p$tmp.datadir, "statistics_loc.bigmemory_description" )
      p$ptr$Stime = file.path( p$tmp.datadir, "statistics_time.bigmemory_description" )

      p$bm =list()
      p$bm$Y =     gsub( "_description", "", basename( p$ptr$Y ) ) 
      p$bm$Ycov =  gsub( "_description", "", basename( p$ptr$Ycov ) ) 
      p$bm$Yloc =  gsub( "_description", "", basename( p$ptr$Yloc ) )
      p$bm$Ytime = gsub( "_description", "", basename( p$ptr$Ytime ) )
      p$bm$P =     gsub( "_description", "", basename( p$ptr$P ) )
      p$bm$Psd =   gsub( "_description", "", basename( p$ptr$Psd ) )
      p$bm$Pn =    gsub( "_description", "", basename( p$ptr$Pn ) )
      p$bm$Pcov =  gsub( "_description", "", basename( p$ptr$Pcov ) )
      p$bm$Ploc =  gsub( "_description", "", basename( p$ptr$Ploc ) )
      p$bm$Ptime = gsub( "_description", "", basename( p$ptr$Ptime ) )
      p$bm$S =     gsub( "_description", "", basename( p$ptr$S ) )
      p$bm$Sloc =  gsub( "_description", "", basename( p$ptr$Sloc ) )
      p$bm$Stime = gsub( "_description", "", basename( p$ptr$Stime ) )
      
      return(p)
    }

    # --------------------------

    if (DS=="save.parameters")  {
      save(p, file=file.path( p$tmp.datadir, "p.rdata") )
      message( "Saved parameters:")
      message( file.path( p$tmp.datadir, "p.rdata") )
    }
    if (DS=="load.parameters")  {
      load(file.path( p$tmp.datadir, "p.rdata") )
      return(p)
    }

    # --------------------------
    if (DS %in% "cleanup" ) {
      p = spacetime_db( p=p, DS="filenames" )
      for (fn in p$ptr ) if (file.exists(fn)) file.remove(fn)
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
      if (exists( "Y_bounds", p) ) {
        bad = which( Y_ < p$Y_bounds[1] | Y_ > p$Y_bounds[2]  )
        if (length( bad) > 0) {
          message( "Y range exceeds that specified in p$Y_bounds ... truncating data") 
          message( p$Y_bounds )
          Y_[ bad] = NaN
        }
      } else {
        p$Y_bounds = range( Y_, na.rm=TRUE) # if not present create it: interpolations must not exceed observed bounds
        message( "Observed Y range:")
        message( p$Y_bounds )
      }

      Y = bigmemory::as.big.matrix( Y_, type="double", 
        backingfile=p$bm$Y, descriptorfile=basename(p$ptr$Y), backingpath=p$tmp.datadir )


     # data coordinates
      Yloc_ = as.matrix( B[, p$variables$LOCS ])
      Yloc = bigmemory::as.big.matrix( Yloc_, type="double", 
        backingfile=p$bm$Yloc, descriptorfile=basename(p$ptr$Yloc), backingpath=p$tmp.datadir )

      # independent variables/ covariate
      if (exists("COV", p$variables)) {
        Ycov_ = as.matrix( B[ , p$variables$COV ] )
        Ycov = bigmemory::as.big.matrix( Ycov_, type="double", 
          backingfile=p$bm$Ycov, descriptorfile=basename(p$ptr$Ycov), backingpath=p$tmp.datadir )
      }

      # prediction times
      if ( exists("TIME", p$variables) ) {
        Ytime_ = as.matrix( B[, p$variables$TIME ] )
        Ytime = bigmemory::as.big.matrix( Ytime_, type="double", 
          backingfile=p$bm$Ytime, descriptorfile=basename(p$ptr$Ytime), backingpath=p$tmp.datadir )
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
        Pcov = bigmemory::as.big.matrix( Pcov_, type="double", 
          backingfile=p$bm$Pcov, descriptorfile=basename(p$ptr$Pcov), backingpath=p$tmp.datadir )
      }
      
      # prediction times
      if (exists("TIME", p$variables)) {
        Ptime_ = as.matrix( B$TIME )
        Ptime = bigmemory::as.big.matrix( Ptime_, type="double", 
          backingfile=p$bm$Ptime, descriptorfile=basename(p$ptr$Ptime), backingpath=p$tmp.datadir )
      }
      
      # predictions and associated stats
      P_ = matrix( NaN, nrow=nrow(B$LOCS), ncol=p$nw*p$ny )
      P = bigmemory::as.big.matrix( P_, type="double", 
        backingfile=p$bm$P, descriptorfile=basename(p$ptr$P), backingpath=p$tmp.datadir )

      Pn = bigmemory::as.big.matrix( Pn_, type="double", 
        backingfile=p$bm$Pn, descriptorfile=basename(p$ptr$Pn), backingpath=p$tmp.datadir )
      
      Psd = bigmemory::as.big.matrix( Psd_, type="double", 
        backingfile=p$bm$Psd, descriptorfile=basename(p$ptr$Psd), backingpath=p$tmp.datadir )

      # prediction coordinates
      Ploc_ = as.matrix( B$LOCS )
      Ploc = bigmemory::as.big.matrix( Ploc_, type="double", 
        backingfile=p$bm$Ploc, descriptorfile=basename(p$ptr$Ploc), backingpath=p$tmp.datadir )

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
      Y = attach.big.matrix( p$ptr$Y )
      hasdata = 1:length(Y)
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA

      # covariates (independent vars)
      if ( exists( "COV", p$variables) ) {
        Ycov = attach.big.matrix( p$ptr$Ycov ) 
        if ( length( p$variables$COV ) == 1 ) {
          bad = which( !is.finite( Ycov[]) )
        } else {
          bad = which( !is.finite( rowSums(Ycov[])) )
        }
        if (length(bad)> 0 ) hasdata[bad] = NA
      }

      ii = na.omit(hasdata)
      ndata = length(ii)
      Yloc = attach.big.matrix( p$ptr$Yloc )
      locs_noise =  runif( ndata*2, min=-p$pres*p$spacetime.noise, max=p$pres*p$spacetime.noise )
      if (!exists("non_convex_hull_alpha", p)) p$non_convex_hull_alpha=20
      boundary=list( polygon = non_convex_hull( Yloc[ii,]+locs_noise, alpha=p$non_convex_hull_alpha, plot=FALSE ) )
      Sloc = attach.big.matrix( p$ptr$Sloc ) # statistical output locations
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
    
    if (DS %in% c( "statistics.initialize", "statistics.size" , "statistics.status", "statistics.box" ) ) {
      p = spacetime_db( p=p, DS="filenames" )  
     
      if ( DS=="statistics.initialize" ) {
        # statistics storage matrix ( aggregation window, coords ) .. no inputs required
        # statistics coordinates
        Sloc_ = as.matrix( expand.grid( p$sbbox$plons, p$sbbox$plats ))
        Sloc = bigmemory::as.big.matrix( Sloc_, type="double", 
          backingfile=p$bm$Sloc, descriptorfile=basename(p$ptr$Sloc), backingpath=p$tmp.datadir )

        S_ = matrix( NaN, nrow=nrow(Sloc_), ncol=length( p$statsvars ) ) # NA forces into logical
        S = bigmemory::as.big.matrix( S_, type="double", 
          backingfile=p$bm$S, descriptorfile=basename(p$ptr$S), backingpath=p$tmp.datadir )
        return( p )
      }

      if ( DS=="statistics.size" ) {
        S = attach.big.matrix( p$ptr$S )
        nS = nrow(S) 
        return(nS)
      }
      
      if (DS == "statistics.box")  {
        sbbox = list( plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$dist.mwin ),
                      plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$dist.mwin ) )
        return(sbbox)
      }


      if ( DS=="statistics.status" ) {
        # find locations for statistic computation and trim area based on availability of data
        # stats:
        S = attach.big.matrix( p$ptr$S )
        bnds = try( spacetime_db( p, DS="boundary" ) )

        if (!is.null(bnds)) {
          if( !("try-error" %in% class(bnds) ) ) {
            # problematic and/or no data (e.g., land, etc.) and skipped
            to.ignore =  which( bnds$inside.polygon == 0 ) # outside boundary
            i = which( is.infinite( S[,1] ) & bnds$inside.polygon != 0 ) # not yet completed (due to a problem)
            j = which( is.nan( S[,1] )  & bnds$inside.polygon != 0 )     # completed
            k = which( is.finite (S[,1])  & bnds$inside.polygon != 0 )   # not yet done
        } } else {
            to.ignore = NA
            i = which( is.infinite( S[,1] )  )  # not yet completed (due to a failed attempt)
            j = which( is.nan( S[,1] )   )      # completed
            k = which( is.finite (S[,1])  )     # not yet done
        }

        out = list(problematic=i, incomplete=j, completed=k, n.total=nrow(S) ,
                     n.incomplete=length(j), n.problematic=length(i), 
                     n.complete=length(k), to.ignore=to.ignore )
        out$prop_incomp=out$n.incomplete / ( out$n.problematic + out$n.incomplete + out$n.complete)
    
        message( paste("Proportion incomplete:", round(out$prop_incomp,5), "\n" )) 
        return( out )
      }
    }


  }

