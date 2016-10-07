
  spacetime_db = function( DS, p, B=NULL, grp=NULL, yr=NULL ) {
    #// usage: low level function to convert data into file-based data obects to permit parallel
    #// data access and manipulation and deletes/updates
    #// B is the xyz or xytz data or the function to get the data to work upon
    #// NaN and Inf are  used instead of NA as ff converts NA initialised matrices as TRue/FALSe
    #// note:: ff seems to treat is.nan and is.na as the same thing 
    #\\ data extraction layer for user objects created by spacetime
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
      # input data stored as a ff file to permit operations with min memory usage
      # split into separate components to reduce filelocking conflicts

      # storage locations for finalized data 
      p$fn = list()
      p$fn$P = file.path( p$savedir, paste( "spacetime", "predictions", "rdata", sep=".") )
      p$fn$S = file.path( p$savedir, paste( "spacetime", "statistics",  "rdata", sep=".") )
      p$fn$results.covar =  file.path( p$project.root, "spacetime", paste( "spatial", "covariance", "rdata", sep=".") )

      p$ptr =list()
      p$ptr$Y =    file.path( p$tmp.datadir, "input.Y.ff")
      p$ptr$Ycov = file.path( p$tmp.datadir, "input.Ycov.ff" )
      p$ptr$Yloc = file.path( p$tmp.datadir, "input.Yloc.ff")
      p$ptr$Ytime = file.path( p$tmp.datadir, "input.Ytime.ff")
      p$ptr$P =    file.path( p$tmp.datadir, "predictions.ff")
      p$ptr$Psd =  file.path( p$tmp.datadir, "predictions_sd.ff")
      p$ptr$Pn =   file.path( p$tmp.datadir, "predictions_n.ff")
      p$ptr$Pcov = file.path( p$tmp.datadir, "predictions_cov.ff")
      p$ptr$Ploc = file.path( p$tmp.datadir, "predictions_loc.ff")
      p$ptr$Ptime = file.path( p$tmp.datadir, "predictions_time.ff")
      p$ptr$S =    file.path( p$tmp.datadir, "statistics.ff")
      p$ptr$Sloc = file.path( p$tmp.datadir, "statistics_loc.ff")
      p$ptr$Stime = file.path( p$tmp.datadir, "statistics_time.ff")
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
      Yb = B[, p$variables$Y ]
      if (exists( "Y_bounds", p) ) {
        bad = which( Yb < p$Y_bounds[1] | Yb > p$Y_bounds[2]  )
        if (length( bad) > 0) {
          message( "Y range exceeds that specified in p$Y_bounds ... truncating data") 
          message( p$Y_bounds )
          Yb[ bad] = NaN
        }
      } else {
        p$Y_bounds = range( Yb, na.rm=TRUE) # if not present create it: interpolations must not exceed observed bounds
        message( "Observed Y range:")
        message( p$Y_bounds )
      }

      Y = ff( Yb, filename=p$ptr$Y, overwrite=TRUE )
      close(Y) # ensure write/sync
      p$ff$Y = Y  # store pointer

      # independent variables/ covariate
      if (exists("COV", p$variables)) {
        Ycovb = as.matrix( B[ , p$variables$COV ] )
        Ycov = ff::ff( Ycovb, dim=dim(Ycovb), filename=p$ptr$Ycov, overwrite=TRUE )
        close(Ycov)
        p$ff$Ycov = Ycov # store pointer
      }

     # data coordinates
      Ylocb = as.matrix( B[, p$variables$LOCS ])
      Yloc = ff::ff( Ylocb, dim=dim(Ylocb), filename=p$ptr$Yloc, overwrite=TRUE )
      close(Yloc)
      p$ff$Yloc = Yloc # store pointer

      # prediction times
      if ( exists("TIME", p$variables) ) {
        Ytimeb = as.matrix( B[, p$variables$TIME ] )
        Ytime = ff::ff( Ytimeb, dim=dim( Ytimeb), filename=p$ptr$Ytime, overwrite=TRUE )
        close(Ytime)
        p$ff$Ytime = Ytime  # store pointer
      }

      return( p ) #return pointers to data
    }


    #---------------------
    if (DS == "predictions.initialize.xy" ) {

      p = spacetime_db( p=p, DS="filenames" ) 
      

      if ( "data.frame" %in% class(B) ) {
        # B = B # nothing to do
      } else if ( "SpatialGridDataFrame" %in% class(B) ) {
        B = slot(B, "data")  
      }
      
    # prediction covariates i.e., independent variables/ covariates
      if (exists("COV", p$variables)) {
        if (is.vector(B$COV) ) {
          Pcovb = as.matrix( B$COV ) 
        } else {
          Pcovb = as.matrix( B$COV[,p$variables$COV ] ) 
        }
        p$ff$Pcov = ff::ff( Pcovb, dim=dim(Pcovb), filename=p$ptr$Pcov, overwrite=TRUE )
        close(p$ff$Pcov)
      }

      # prediction coordinates
      Plocb = as.matrix( B[, p$variables$LOCS ] )
      Ploc = ff::ff( Plocb, dim=dim(Plocb), filename=p$ptr$Ploc, overwrite=TRUE )
      close(Ploc)
      p$ff$Ploc = Ploc  # store pointer
  
      # predictions and associated stats
      Pb = matrix( NaN, nrow=p$nPreds, ncol=3 ) # pred, count, sd
      P = ff::ff( Pb, dim=dim(P), filename=p$ptr$P, overwrite=TRUE )
      close(P)
      p$ff$P = P  # store pointer
      
      return( p )
    }



    #---------------------
    if (DS == "predictions.initialize.xyt" ) {
      # prediction covariates i.e., independent variables/ covariates

      p = spacetime_db( p=p, DS="filenames" ) 
  
      if (exists("COV", p$variables)) {
        if (is.vector(B$COV) ) {
          Pcovb = as.matrix( B$COV ) 
        } else {
          Pcovb = as.matrix( B$COV[,p$variables$COV ] ) 
        }
        p$ff$Pcov = ff::ff( Pcovb, dim=dim(Pcovb), filename=p$ptr$Pcov, overwrite=TRUE )
        close(p$ff$Pcov)
      }
          
      # predictions and associated stats
      # predictions and associated stats
      P = ff::ff( NaN, dim=c( nrow(B$LOCS), p$ny), filename=p$ptr$P, overwrite=TRUE )
      close(P)      
      p$ff$P = P  # store pointer

      Pn = ff::ff( Pb, dim=dim(Pb), filename=p$ptr$Pn, overwrite=TRUE )
      close(Pn)      
      p$ff$Pn = Pn  # store pointer
      
      Psd = ff::ff( Pb, dim=dim(Pb), filename=p$ptr$Psd, overwrite=TRUE )
      close(Psd)      
      p$ff$Psd = Psd  # store pointer

      # prediction coordinates
      Plocb = as.matrix( B[, p$variables$LOCS ] )
      Ploc = ff::ff( Plocb, dim=dim(Plocb), filename=p$ptr$Ploc, overwrite=TRUE )
      close(Ploc)
      p$ff$Ploc = Ploc  # store pointer

      # prediction times
      Ptime_ = as.matrix( B[, p$variables$TIME ] )
      Ptime = ff::ff( Ptime_, dim=dim(Ptime_), filename=p$ptr$Ptime, overwrite=TRUE )
      close(Ptime)
      p$ff$Ptime = Ptime  # store pointer
         
      return( p )
    }


    #---------------------
    if (DS == "predictions.initialize.xyts" ) {

      p = spacetime_db( p=p, DS="filenames" ) 
      
      if (exists("COV", p$variables)) {
        if (is.vector(B$COV) ) {
          Pcovb = as.matrix( B$COV ) 
        } else {
          Pcovb = as.matrix( B$COV[,p$variables$COV ] ) 
        }
        p$ff$Pcov = ff::ff( Pcovb, dim=dim(Pcovb), filename=p$ptr$Pcov, overwrite=TRUE )
        close(p$ff$Pcov)
      }
      
      # predictions and associated stats
      P = ff::ff( NaN, dim=c( nrow(B$LOCS), p$nw*p$ny), filename=p$ptr$P, overwrite=TRUE )
      close(P)      
      p$ff$P = P  # store pointer

      Pn = ff::ff(NaN, dim=dim(P), filename=p$ptr$Pn, overwrite=TRUE )
      close(Pn)      
      p$ff$Pn = Pn  # store pointer
      
      Psd = ff::ff( NaN, dim=dim(P), filename=p$ptr$Psd, overwrite=TRUE )
      close(Psd)      
      p$ff$Psd = Psd  # store pointer

      # prediction coordinates
      Plocb = as.matrix( B$LOCS )
      Ploc = ff::ff( Plocb, dim=dim(Plocb), filename=p$ptr$Ploc, overwrite=TRUE )
      close(Ploc)
      p$ff$Ploc = Ploc  # store pointer

      # prediction times
      Ptime_ = as.matrix( B$TIME )
      Ptime = ff::ff( Ptime_, dim=dim(Ptime_), filename=p$ptr$Ptime, overwrite=TRUE )
      close(Ptime)
      p$ff$Ptime = Ptime  # store pointer

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
      Y = p$ff$Y
      hasdata = 1:length(Y)
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA
      close(Y)

      # covariates (independent vars)
      if ( exists( "COV", p$variables) ) {
        Ycov = p$ff$Ycov 
        if ( length( p$variables$COV ) == 1 ) {
          bad = which( !is.finite( Ycov[]) )
        } else {
          bad = which( !is.finite( rowSums(Ycov[])) )
        }
        if (length(bad)> 0 ) hasdata[bad] = NA
        close(Ycov)
      }

      ii = na.omit(hasdata)
      ndata = length(ii)

      Yloc = p$ff$Yloc 
      locs_noise =  runif( ndata*2, min=-p$pres*p$spacetime.noise, max=p$pres*p$spacetime.noise )
      boundary=list( polygon = non_convex_hull( Yloc[ii,]+locs_noise, alpha=p$non_convex_hull_alpha, plot=FALSE ) )
      Sloc = p$ff$Sloc  # statistical output locations
      boundary$inside.polygon = point.in.polygon( Sloc[,1], Sloc[,2],
          boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )
      
      save( boundary, file=fn, compress=TRUE )
      plot( Yloc[ii,], pch=".", col="grey" ) # data locations
      points( Sloc[which(boundary$inside.polygon==1),], pch=".", col="orange" )
      lines( boundary$polygon[] , col="green", pch=2 )
      message( "Check the map of data and boundaries. ")
      message( "If not suitable, set another value for p$non_convex_hull_alpha value (radius; distance) ")
      message( "and re-run spacetime() with the flag: overwrite=FALSE " )
      close(Sloc)
      close(Yloc)      
      return( fn )
    }


    # -----------------
    
    if (DS %in% c( "statistics.initialize", "statistics.size" , "statistics.status", "statistics.box" ) ) {
      p = spacetime_db( p=p, DS="filenames" )  
     
      if ( DS=="statistics.initialize" ) {
        # statistics storage matrix ( aggregation window, coords ) .. no inputs required
        # statistics coordinates
        Sloc_ = as.matrix( expand.grid( p$sbbox$plons, p$sbbox$plats ))
        Sloc = ff::ff( Sloc_, dim=dim(Sloc_), filename=p$ptr$Sloc, overwrite=TRUE )
        close(Sloc)
        p$ff$Sloc = Sloc

        S_ = matrix( NaN, nrow=nrow(Sloc_), ncol=length( p$statsvars ) ) # NA forces into logical
        S = ff::ff( S_, dim=dim(S_), filename=p$ptr$S, overwrite=TRUE )
        close(S)
        p$ff$S = S

        return( p )
      }
      

      if ( DS=="statistics.size" ) {
        nS = nrow(p$ff$S) 
        close(S) 
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
        S = p$ff$S
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
        close(S)
        cat( paste("Proportion incomplete:", round(out$prop_incomp,5), "\n" )) 
        return( out )
      }
    }


  }

