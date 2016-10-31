
  spacetime_db = function( DS, p, B=NULL, yr=NULL, ret="mean"  ) {
    #// usage: low level function to convert data into file-based data obects to permit parallel
    #// data access and manipulation and deletes/updates
    #// B is the xyz or xytz data or the function to get the data to work upon

    # --------------------------

    if (DS %in% "filenames" ) {
      # input data stored as a bigmemory file to permit operations with min memory usage
      # split into separate components to reduce filelocking conflicts
         
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
      p$cache$Ylogit = file.path( p$stloc, "Ylogit.cache" )

      p$cache$P0 =    file.path( p$stloc, "predictions0.cache" ) # offsets from covar model
      p$cache$P0sd =  file.path( p$stloc, "predictions0sd.cache" ) # offsets from covar model
      p$cache$P =     file.path( p$stloc, "predictions.cache" )
      p$cache$Psd =   file.path( p$stloc, "predictions_sd.cache" )
      p$cache$Pn =    file.path( p$stloc, "predictions_n.cache" )
      p$cache$Pcov =  file.path( p$stloc, "predictions_cov.cache" )
      p$cache$Ploc =  file.path( p$stloc, "predictions_loc.cache" )
      p$cache$Ptime = file.path( p$stloc, "predictions_time.cache" )
      p$cache$Plogit = file.path( p$stloc, "Plogit.cache" )
      p$cache$Plogitsd = file.path( p$stloc, "Plogitsd.cache" )

      p$cache$S =     file.path( p$stloc, "statistics.cache" )
      p$cache$Sloc =  file.path( p$stloc, "statistics_loc.cache" )
      p$cache$Stime = file.path( p$stloc, "statistics_time.cache" )
      p$cache$Sflag =     file.path( p$stloc, "statistics_flag.cache" )

      p$cache$Mat2Ploc = file.path( p$stloc, "Mat2Ploc.cache" )

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
      load( file.path( p$stloc, "p.rdata") )
      return(p)
    }

    # --------------------------
    if (DS %in% "cleanup" ) {
      for (fn in p$cache ) if (file.exists(fn)) file.remove(fn)
      for (fn in p$bm ) if (file.exists(fn)) file.remove(fn)
      return( "done" )
    }

    # -----------------
    
    if (DS %in% c( "statistics.status", "statistics.box" ) ) {
          
      if (DS == "statistics.box")  {
        sbbox = list( plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$spacetime_distance_statsgrid ),
                      plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$spacetime_distance_statsgrid ) )
        return(sbbox)
      }

      if ( DS=="statistics.status" ) {
        # find locations for statistic computation and trim area based on availability of data
        # stats:
        Sflag = spacetime_attach( p$storage.backend, p$ptr$Sflag )
        
        i = which( is.infinite( Sflag[] )  )  # not yet completed (due to a failed attempt)
        j = which( is.nan( Sflag[] )   )      # incomplete
        k = which( is.finite (Sflag[] )  )     # completed
        bnds = try( spacetime_db( p, DS="boundary" ) )
        if (!is.null(bnds)) {
          if( !("try-error" %in% class(bnds) ) ) {
            if (0) {
              # to reset the flags
              i = which( is.infinite( Sflag[] )  )  # not yet completed (due to a failed attempt)
              Sflag[i] = NaN
  
              l =  which( bnds$inside.polygon == 0 ) # outside boundary
              if (length(l)>0) Sflag[l] = Inf  # outside of data area
            }
        }}

        out = list(problematic=i, incomplete=j, completed=k, n.total=length(Sflag) ,
                     n.incomplete=length(j), n.problematic=length(i), 
                     n.complete=length(k) )
        out$prop_incomp=out$n.incomplete / ( out$n.incomplete + out$n.complete)
        message( paste("Proportion incomplete:", round(out$prop_incomp,5), "\n" )) 
        return( out )
      }
    }

   
    #---------------------
   
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
  
    if (DS %in% c("model.covariates", "model.covariates.redo") ) {
      
      fn.covmodel =  file.path( p$project.root, "spacetime", paste( "spatial", "covariate.model", p$spacetime_engine, p$spacetime_family, "rdata", sep=".") )

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
          theta=p$spacetime_distance_statsgrid, xwidth=p$spacetime_distance_statsgrid*10, ywidth=p$spacetime_distance_statsgrid*10 )
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

    #-------------

    if (DS=="presence.absense") {

      Y = spacetime_attach( p$storage.backend, p$ptr$Y )
      z = which( Y == 0) # assumed to be real zeros
      i = which( Y >  0)  # positive values
      
      # determine quantiles
      Yq = rep( 0, length(Y) )
      Yq[z] = 1
    
      pr = ecdf(Y[i])( Y[i] )
      ix = which( pr ==1 )
      if ( !( length(ix) %in% c(0, length(x)) ))  pr[ix] = max( pr[-ix] )
      Yq[i] = pr

      s01 = which( Yq < p$habitat.threshold.quantile )  # buffer zone
      s0 = unique( c(s01, sz ) )
      s1 = which( Yq >= p$habitat.threshold.quantile )

      # determine presence-absence
      YPA =  rep( NA, length(Y) )
      YPA[s1] = 1  
      YPA[s0] = 0
      YPA[z] = 0

      return(YPA)
    }

  }


