spacetime_model_predict = function( p, Si, YiU, pa ) {

  res =NULL
  Yloc = p$ptr$Yloc

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
      Sloc = p$ptr$Sloc
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
            res$spacetime_stats["ar_1"] = coef( lm( pac$mean[1:(length(piid) - 1)] ~ pac$mean[2:(length(piid))] + 
0 ) )
          }
        }          
      }
    }
    return(res)
  }

  if (p$spacetime_engine=="kernel.density") res = spacetime__kerneldensity(  )
  if (p$spacetime_engine=="LaplacesDemon") res = spacetime__LaplacesDemon(  )
  if (p$spacetime_engine=="inla") res = spacetime__inla()


}
