
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

      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$bm = p$cache
        for (i in 1:length(p$bm)) p$bm[[i]] = gsub(".cache$", ".bigmemory", p$bm[[i]] )
      }

      if (p$storage.backend == "bigmemory.ram" ) {
        p$bm=list() # initial storage of ram objects
      }      
      
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
      for (fn in p$cache ) if (length(fn)>0) if (file.exists(fn)) file.remove(fn)
      for (fn in p$bm ) if (length(fn)>0)  if (file.exists(fn)) file.remove(fn)
      return( "done" )
    }

    # -----------------
    
    if (DS %in% c( "statistics.status", "statistics.box", "statistics.reset.problem.locations" ) ) {
          
      if (DS == "statistics.box")  {
        sbbox = list( plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$spacetime_distance_statsgrid ),
                      plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$spacetime_distance_statsgrid ) )
        return(sbbox)
      }

      if ( DS=="statistics.status" ) {
        # find locations for statistic computation and trim area based on availability of data
        # stats:
        bnds = try( spacetime_db( p, DS="boundary" ) )
        ioutside = NA
        if (!is.null(bnds)) {
          if( !("try-error" %in% class(bnds) ) ) {
            ioutside = which( bnds$inside.polygon == 0 ) # outside boundary
        }}
        Sflag = spacetime_attach( p$storage.backend, p$ptr$Sflag )
        itodo = setdiff( which( is.nan( Sflag[] )), ioutside)       # incomplete
        idone = setdiff( which( is.finite (Sflag[] )  ), ioutside)      # completed
        iskipped = which( is.infinite( Sflag[] )  ) # skipped due to problems or out of bounds
        iproblems = setdiff( iskipped, ioutside)    # not completed due to a failed attempt
        out = list(problematic=iproblems, skipped=iskipped, todo=itodo, completed=idone, outside=ioutside,
                   n.total=length(Sflag) , n.skipped=length(iskipped),
                   n.todo=length(itodo), n.problematic=length(iproblems), 
                   n.outside=length(which(is.finite(ioutside))),
                   n.complete=length(idone) )
        out$prop_incomp=out$n.todo / ( out$n.todo + out$n.complete)
        message( paste("Proportion to do:", round(out$prop_incomp,5), "\n" )) 
        return( out )
      }
      
      if ( DS=="statistics.reset.problem.locations" ) {
        # to reset all rejected locations 
        Sflag = spacetime_attach( p$storage.backend, p$ptr$Sflag )
        o = spacetime_db( p, DS="statistics.status" )
        if (length(which(is.finite(o$skipped))) > 0) Sflag[o$skipped] = NaN  # to reset all the flags
        if (length(which(is.finite(o$outside))) > 0) Sflag[o$outside] = Inf  # flag area outside of data boundary to skip
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
      
      fn.covmodel =  file.path( p$project.root, "spacetime", paste( "spatial", "covariate.model", p$spacetime_covariate_modeltype, "rdata", sep=".") )

      if (DS =="model.covariates") {
        covmodel = NULL
        if (file.exists( fn.covmodel ))  load(fn.covmodel)
        return(covmodel)
      }  

       good = which( is.finite (rowSums(B[ , c(p$variables$Y,p$variables$COV) ])) )
   
      # as a first pass, model the time-independent factors as a user-defined model
      if (p$spacetime_covariate_modeltype=="gam") {
        covmodel = try( 
          gam( p$spacetime_covariate_modelformula, data=B, optimizer=c("outer","bfgs"), family=p$spacetime_family ) ) 
      }

      if (p$spacetime_covariate_modeltype=="bayesx") {
        if ( !exists( "bayesx_covariate_method", p) ) p$bayesx_covariate_method="REML"  # slightly more smoothing than the REML method
        covmodel = try( 
          bayesx( p$spacetime_covariate_modelformula, data=B, family=p$bayesx_covariate_family,  method=p$bayesx_covariate_method, na.action="na.omit" ) ) 
      }

      if ( "try-error" %in% class(covmodel) ) stop( "The covariate model was problematic" )
      print( summary( covmodel ) )
      save( covmodel, file= fn.covmodel, compress=TRUE )
    
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
      

      if (exists("model.covariates.globally", p) && p$model.covariates.globally ) {
        P0 = spacetime_attach( p$storage.backend, p$ptr$P0 )
        P0sd = spacetime_attach( p$storage.backend, p$ptr$P0sd )
        PP = PP + P0 
        PPsd = sqrt( P0sd^2 + PPsd^2) # simpleadditive independent errors assumed
      }

      if ( exists("TIME", p)) {
        # outputs are on yearly breakdown
        for ( r in 1:length(p$spacetime_yrs) ) {
          y = p$spacetime_yrs[r]
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
          P = PP
          V = PPsd
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
    
      Ploc = spacetime_attach( p$storage.backend, p$ptr$Ploc )
  
      S = spacetime_attach( p$storage.backend, p$ptr$S )
      Sloc = spacetime_attach( p$storage.backend, p$ptr$Sloc )
    
      # locations of the new (output) coord system
      locsout = expand.grid( p$plons, p$plats ) # final output grid
      attr( locsout , "out.attrs") = NULL
      names( locsout ) = p$variables$LOCS

      stats = matrix( NaN, ncol=length( p$statsvars ), nrow=nrow( locsout) )  # output data
      colnames(stats)=p$statsvars

      # map of row, col indices of input data in the new (output) coordinate system
      l2M = cbind( ( Sloc[,1]-p$plons[1])/p$pres + 1, (Sloc[,2]-p$plats[1])/p$pres + 1)
     
      # matrix representation of the output surface
      M = matrix( NA, nrow=p$nplons, ncol=p$nplats) 
      
      for ( i in 1:length( p$statsvars ) ) {
        M = M[] * NA  # init
        M[l2M] = S[,i] # fill with data in correct locations
        # Z = fields::interp.surface( list( x=p$plons, y=p$plats, z=M ), loc=locsout )
        # ii = which( !is.finite( Z ) )
        # if ( length( ii) > 0 ) {
          # try again ..
          # Z[ii] = fields::interp.surface( list( x=p$plons, y=p$plats, z=M ), loc=locsout[ii,] )
        # }
        # ii = which( !is.finite( Z ) )
        #if ( length( ii) > 0 ) {
        Z =  fields::image.smooth( M, dx=p$pres, dy=p$pres, theta=p$theta )$z  
        Z[l2M] = S[,i]  # return raw data to predicted surface
        #  Z[ii] = Zii[ii]
        # }
        stats[,i] = Z
      }

      boundary = try( spacetime_db( p, DS="boundary" ) )
      if (!is.null(boundary)) {
        if( !("try-error" %in% class(boundary) ) ) {
        inside.polygon = point.in.polygon( locsout[,1], locsout[,2],
          boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )
        o = which( inside.polygon == 0 ) # outside boundary
        if (length(o) > 0) stats[o,] = NA         
      }}

      # subset to match to Ploc
      Ploc_id = paste( Ploc[,1], Ploc[,2], sep="~" )
      locsout_id = paste( locsout$plon, locsout$plat, sep="~" )
      good = match( Ploc_id, locsout_id )

      bad = which( !is.finite(good))
      if (length(bad) > 0 ) good = good[-bad]
      stats = stats[ good, ]

      save( stats, file=p$fn.S, compress=TRUE )

      # lattice::levelplot( Pstats[,1] ~ Ploc[,1]+Ploc[,2])
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
      Ybin =  rep( NA, length(Y) )
      Ybin[s1] = 1  
      Ybin[s0] = 0
      Ybin[z] = 0

      return(Ybin)
    }

  }


