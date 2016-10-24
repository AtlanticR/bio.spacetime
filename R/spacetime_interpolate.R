
spacetime_interpolate = function( ip=NULL, p ) {
  #\\ core function to intepolate (model and predict) in parllel

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  #---------------------
  # data for modelling
  # dependent vars # already link-transformed in spacetime_db("dependent")
  S = spacetime_attach( p$storage.backend, p$ptr$S )
  Sflag = spacetime_attach( p$storage.backend, p$ptr$Sflag )
  
  Sloc = spacetime_attach( p$storage.backend, p$ptr$Sloc )
  Yloc = spacetime_attach( p$storage.backend, p$ptr$Yloc )
  Y = spacetime_attach( p$storage.backend, p$ptr$Y )
  Yloc = spacetime_attach( p$storage.backend, p$ptr$Yloc )

  if ( p$storage.backend != "bigmemory.ram" ) {
    # force copy into RAM
    Sloc = Sloc[]
    Yloc = Yloc[]
    Y = Y[]
    Yloc = Yloc[]
  }

  Yi = spacetime_attach( p$storage.backend, p$ptr$Yi )
  Yi = as.vector(Yi[])  #force copy to RAM

  # main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    Si = p$runs[ iip, "locs" ]
    if ( is.infinite( Sflag[Si] ) ) next() 
    if ( !is.nan( Sflag[Si] ) ) next() 
    Sflag[Si] = Inf   # over-written below if successful else if a run fails it does not get revisited 
    print( iip )

    # find data withing a given distance / number 
    pib = point_in_block( Sloc=Sloc, Si=Si, Yloc=Yloc, Yi=Yi, 
      dist.max=p$dist.max, dist.min=p$dist.min, n.min=p$n.min, n.max=p$n.max, 
      upsampling=p$upsampling, downsampling=p$downsampling, resize=TRUE ) 
    
    if ( is.null(pib)) next()
    
    dist.cur = pib$dist
    ndata = length(pib$U)
    if ((ndata < p$n.min) | (ndata > p$n.max) ) next()
    YiU = Yi[pib$U]  
    rm(pib)
    
    # construct prediction/output grid area ('pa')
    pa = NULL
    pa = spacetime_prediction_area( p, Si, dist.cur ) 
    if (is.null(pa)) next()
    
    # model and prediction
    res =NULL
  
    if ( p$spacetime_engine %in% 
      c( "harmonics.1", "harmonics.2", "harmonics.3", "harmonics.1.depth",
         "seasonal.basic", "seasonal.smoothed", "annual", "gam"  ) ) {
      res = NULL
      res = spacetime__harmonics( p, YiU, Si, pa )
      if ( is.null(res)) return(NULL)     

      if (exists("spacetime_variogram_engine", p) ) {
        sp.stat = NULL
        sp.stat = try( spacetime_variogram(  Yloc[YiU,], Y[YiU], methods=p$spacetime_variogram_engine ) )
        if (!is.null(sp.stat) && !("try-error" %in% class(sp.stat)) ){
          res$spacetime_stats["sdSpatial"] = sqrt( sp.stat[[p$spacetime_variogram_engine]]$varSpatial )
          res$spacetime_stats["sdObs"] = sqrt( sp.stat[[p$spacetime_variogram_engine]]$varObs )
          res$spacetime_stats["range"] = sp.stat[[p$spacetime_variogram_engine]]$range
          res$spacetime_stats["phi"] = sp.stat[[p$spacetime_variogram_engine]]$phi
          res$spacetime_stats["nu"] = sp.stat[[p$spacetime_variogram_engine]]$nu
        }
      }

      if ( exists("TIME", p$variables) ){
        # annual ts, seasonally centered and spatially 
        # pa_i = which( Sloc[Si,1]==Ploc[,1] & Sloc[Si,2]==Ploc[,2] )
        Sloc = spacetime_attach( p$storage.backend, p$ptr$Sloc )
        pac_i = which( res$predictions$plon==Sloc[Si,1] & res$predictions$plat==Sloc[Si,2] )
        if (length(pac_i) > 5) {
          pac = res$predictions[ pac_i, ]
          pac$dyr = pac[, p$variables$TIME] - trunc(pac[, p$variables$TIME] )
          piid = which( zapsmall( pac$dyr - p$dyear_centre) == 0 )
          pac = pac[ piid, c(p$variables$TIME, "mean")]
          pac = pac[ order(pac[,p$variables$TIME]),]
          if (length(piid) > 5 ) {
            ts.stat = NULL
            ts.stat = try( spacetime_timeseries( pac$mean, method="fft" ) )
            if (!is.null(ts.stat) && !("try-error" %in% class(ts.stat)) ){
              res$spacetime_stats["ar_timerange"] = ts.stat$quantilePeriod 
              if (length(which (is.finite(pac$mean))) > 5 ) {
                ar1 = try( lm( pac$mean[1:(length(piid) - 1)] ~ pac$mean[2:(length(piid))] + 0, na.action="na.omit") )
                res$spacetime_stats["ar_1"] = coef( ar1 )
              }
            }
          }          
        }
      }
      return(res)
    }

    if (p$spacetime_engine=="kernel.density") res = spacetime__kerneldensity(  )
    if (p$spacetime_engine=="LaplacesDemon") res = spacetime__LaplacesDemon(  )
    if (p$spacetime_engine=="inla") res = spacetime__inla()



    if ( is.null(res)) next()
    rm(pa)
    
    for ( k in 1: length(p$statsvars) ) {
      if (exists( p$statsvars[k], res$spacetime_stats )) {
        S[Si,k] = res$spacetime_stats[[ p$statsvars[k] ]]
      }
    }

    spacetime_predictions_save( p, res$predictions ) # update P (predictions) .. slow!! .. try diff cache size for P, Pn Psd

      if (0) {
        v = res$predictions
        if ( exists("TIME", p$variables) ){
          v = v[which( v[,p$variables$TIME]==1990.55),]
        }
        require(lattice)
        levelplot( mean ~ plon+plat, v, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
      }

   
    # ----------------------
    # save statistics: do last. it is an indicator of completion of all tasks 
    # .. restarts would be broken otherwise
    Sflag[Si] = 1  # done .. any finite value

  }  # end for loop

}



