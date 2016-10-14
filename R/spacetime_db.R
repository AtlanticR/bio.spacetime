
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
      
      # p$ptr =list()
      # p$ptr$Y0 =    file.path( p$tmp.datadir, "input.Y0.bigmemory_description" ) # raw data
      # p$ptr$Y =     file.path( p$tmp.datadir, "input.Y.bigmemory_description" ) # residuals of covar model or raw data if none
      # p$ptr$Ycov =  file.path( p$tmp.datadir, "input.Ycov.bigmemory_description"  )
      # p$ptr$Yloc =  file.path( p$tmp.datadir, "input.Yloc.bigmemory_description" )
      # p$ptr$Ytime = file.path( p$tmp.datadir, "input.Ytime.bigmemory_description" )
      # p$ptr$P0 =    file.path( p$tmp.datadir, "predictions0.bigmemory_description" ) # offsets from covar model
      # p$ptr$P =     file.path( p$tmp.datadir, "predictions.bigmemory_description" )
      # p$ptr$Psd =   file.path( p$tmp.datadir, "predictions_sd.bigmemory_description" )
      # p$ptr$Pn =    file.path( p$tmp.datadir, "predictions_n.bigmemory_description" )
      # p$ptr$Pcov =  file.path( p$tmp.datadir, "predictions_cov.bigmemory_description" )
      # p$ptr$Ploc =  file.path( p$tmp.datadir, "predictions_loc.bigmemory_description" )
      # p$ptr$Ptime = file.path( p$tmp.datadir, "predictions_time.bigmemory_description" )
      # p$ptr$S =     file.path( p$tmp.datadir, "statistics.bigmemory_description" )
      # p$ptr$Sloc =  file.path( p$tmp.datadir, "statistics_loc.bigmemory_description" )
      # p$ptr$Stime = file.path( p$tmp.datadir, "statistics_time.bigmemory_description" )


      p$ptr_data =list()
      p$ptr_data$Y0 =    gsub( "_description", "",  p$ptr$Y0 ) 
      p$ptr_data$Y =     gsub( "_description", "",  p$ptr$Y ) 
      p$ptr_data$Ycov =  gsub( "_description", "",  p$ptr$Ycov ) 
      p$ptr_data$Yloc =  gsub( "_description", "",  p$ptr$Yloc )
      p$ptr_data$Ytime = gsub( "_description", "",  p$ptr$Ytime )
      p$ptr_data$P0 =    gsub( "_description", "",  p$ptr$P0 )
      p$ptr_data$P =     gsub( "_description", "",  p$ptr$P )
      p$ptr_data$Psd =   gsub( "_description", "",  p$ptr$Psd )
      p$ptr_data$Pn =    gsub( "_description", "",  p$ptr$Pn )
      p$ptr_data$Pcov =  gsub( "_description", "",  p$ptr$Pcov )
      p$ptr_data$Ploc =  gsub( "_description", "",  p$ptr$Ploc )
      p$ptr_data$Ptime = gsub( "_description", "",  p$ptr$Ptime )
      p$ptr_data$S =     gsub( "_description", "",  p$ptr$S )
      p$ptr_data$Sloc =  gsub( "_description", "",  p$ptr$Sloc )
      p$ptr_data$Stime = gsub( "_description", "",  p$ptr$Stime )

      p$bm =list()
      p$bm$Y0 =    basename( p$ptr_data$Y0 ) 
      p$bm$Y =     basename( p$ptr_data$Y ) 
      p$bm$Ycov =  basename( p$ptr_data$Ycov ) 
      p$bm$Yloc =  basename( p$ptr_data$Yloc )
      p$bm$Ytime = basename( p$ptr_data$Ytime )
      p$bm$P0 =    basename( p$ptr_data$P0 )
      p$bm$P =     basename( p$ptr_data$P )
      p$bm$Psd =   basename( p$ptr_data$Psd )
      p$bm$Pn =    basename( p$ptr_data$Pn )
      p$bm$Pcov =  basename( p$ptr_data$Pcov )
      p$bm$Ploc =  basename( p$ptr_data$Ploc )
      p$bm$Ptime = basename( p$ptr_data$Ptime )
      p$bm$S =     basename( p$ptr_data$S )
      p$bm$Sloc =  basename( p$ptr_data$Sloc )
      p$bm$Stime = basename( p$ptr_data$Stime )
      
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
      for (fn in p$ptr_data ) if (file.exists(fn)) file.remove(fn)
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
      p$ptr$Y0 = ff( Y_, dim=dim(Y_), file=p$ptr_data$Y0, overwrite=TRUE )

      if (exists( "spacetime_covariate_spacetime_engine_modelformula", p)) {
        Yresid_ = residuals( spacetime_db( p=p, DS="model.covariates") )
      } else {
        Yresid_ = Y_  # if no residual model simply use the raw data
      }
      p$ptr$Y = ff( Yresid_, dim=dim(Yresid_), file=p$ptr_data$Y, overwrite=TRUE )

     # data coordinates
      Yloc_ = as.matrix( B[, p$variables$LOCS ])
      p$ptr$Yloc = ff( Yloc_, dim=dim(Yloc_), file=p$ptr_data$Yloc, overwrite=TRUE )

      # independent variables/ covariate
      if (exists("COV", p$variables)) {
        Ycov_ = as.matrix( B[ , p$variables$COV ] )
        p$ptr$Ycov = ff( Ycov_, dim=dim(Ycov_), file=p$ptr_data$Ycov, overwrite=TRUE )
      }

      # data times
      if ( exists("TIME", p$variables) ) {
        Ytime_ = as.matrix( B[, p$variables$TIME ] )
        p$ptr$Ytime = ff( Ytime_, dim=dim(Ytime_), file=p$ptr_data$Ytime, overwrite=TRUE )
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
        p$ptr$Pcov = ff( Pcov_, dim=dim(Pcov_), file=p$ptr_data$Pcov, overwrite=TRUE )
      }
      
      # prediction times
      if (exists("TIME", p$variables)) {
        Ptime_ = as.matrix( B$TIME )
        p$ptr$Ptime = ff( Ptime_, dim=dim(Ptime_), file=p$ptr_data$Ptime, overwrite=TRUE )
      }
      
      # predictions and associated stats
      P_ = matrix( NaN, nrow=nrow(B$LOCS), ncol=p$nw*p$ny )
      p$ptr$P = ff( P_, dim=dim(P_), file=p$ptr_data$P, overwrite=TRUE )
      p$ptr$Pn = ff( P_, dim=dim(P_), file=p$ptr_data$Pn, overwrite=TRUE )
      p$ptr$Psd = ff( P_, dim=dim(P_), file=p$ptr_data$Psd, overwrite=TRUE )

      # prediction coordinates
      Ploc_ = as.matrix( B$LOCS )
      p$ptr$Ploc = ff( Ploc_, dim=dim(Ploc_), file=p$ptr_data$Ploc, overwrite=TRUE )

      rm(P_);gc()

      # prediction offsets from an additive (global) model of covariates (if any, zero valued otherwise) 
      if ( exists("COV", p$variables) && exists("spacetime_covariate_spacetime_engine_modelformula", p ) ) {
        # assuming model is correct..
        covmodel = spacetime_db( p=p, DS="model.covariates") 
        covars = data.frame(B$COV)
        names(covars) = p$variables$COV
        Pcov = predict( covmodel, newdata=covars, type="response", se.fit=T ) 
        rm (covmodel, covars);gc()
        p$ptr$P0_ = ff( Pcov$fit, dim=dim(Pcov$fit), file=p$ptr_data$P0, overwrite=TRUE )
        p$ptr$P0sd_ = ff( Pcov$se.fit, dim=dim(Pcov$se.fit), file=p$ptr_data$P0sd, overwrite=TRUE )
      } else {
        P0_ = matrix( 0, nrow=nrow(B$LOCS), ncol=p$nw*p$ny )
        p$ptr$P0 = ff( P0_, dim=dim(P0_), file=p$ptr_data$P0, overwrite=TRUE )
        p$ptr$P0sd = ff( P0sd_, dim=dim(P0sd_), file=p$ptr_data$P0sd, overwrite=TRUE )
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
      Y = p$ptr$Y 
      hasdata = 1:length(Y)
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA

      # covariates (independent vars)
      if ( exists( "COV", p$variables) ) {
        Ycov = p$ptr$Ycov 
        if ( length( p$variables$COV ) == 1 ) {
          bad = which( !is.finite( Ycov[]) )
        } else {
          bad = which( !is.finite( rowSums(Ycov[])) )
        }
        if (length(bad)> 0 ) hasdata[bad] = NA
      }

      ii = na.omit(hasdata)
      ndata = length(ii)
      Yloc = p$ptr$Yloc
      locs_noise =  runif( ndata*2, min=-p$pres*p$spacetime.noise, max=p$pres*p$spacetime.noise )
      if (!exists("non_convex_hull_alpha", p)) p$non_convex_hull_alpha=20
      boundary=list( polygon = non_convex_hull( Yloc[ii,]+locs_noise, alpha=p$non_convex_hull_alpha, plot=FALSE ) )
      Sloc = p$ptr$Sloc # statistical output locations
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
        p$ptr$Sloc = ff( Sloc_, dim=dim(Sloc_), filename=p$ptr_data$Sloc, overwrite=TRUE )

        S_ = matrix( NaN, nrow=nrow(Sloc_), ncol=length( p$statsvars ) ) # NA forces into logical
        p$ptr$S = ff( S_, dim=dim(S_), file=p$ptr_data$S, overwrite=TRUE )
        return( p )
      }

      if ( DS=="statistics.size" ) {
        S = p$ptr$S 
        nS = nrow(S) 
        return(nS)
      }
      
      if (DS == "statistics.box")  {
        sbbox = list( plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$spacetime.prediction.dist.min ),
                      plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$spacetime.prediction.dist.min ) )
        return(sbbox)
      }

      if ( DS=="statistics.status" ) {
        # find locations for statistic computation and trim area based on availability of data
        # stats:
        S = p$ptr$S 
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

      PP =  p$ptr$P 
      PPsd =  p$ptr$Psd 
      P0 =  p$ptr$P0 
      P0sd =  p$ptr$P0sd 

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

      S =  p$ptr$S 
      Sloc =  p$ptr$Sloc 
    
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
      Ploc =  p$ptr$Ploc
      
      bad = which( !is.finite(pa$i))
      if (length(bad) > 0 ) pa = pa[-bad,]
      if (nrow(pa)< 5) next()
      pa$plon = Ploc[ pa$i, 1]
      pa$plat = Ploc[ pa$i, 2]

      save( stats,  file=p$fn.S, compress=TRUE )

    }

  }
