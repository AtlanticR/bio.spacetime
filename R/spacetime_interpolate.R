
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
  Ploc = spacetime_attach( p$storage.backend, p$ptr$Ploc )
  Yloc = spacetime_attach( p$storage.backend, p$ptr$Yloc )

  Y = spacetime_attach( p$storage.backend, p$ptr$Y )

  P = spacetime_attach( p$storage.backend, p$ptr$P )
  Pn = spacetime_attach( p$storage.backend, p$ptr$Pn )
  Psd = spacetime_attach( p$storage.backend, p$ptr$Psd )

  if (exists("COV", p$variables)) {
    Ycov = spacetime_attach( p$storage.backend, p$ptr$Ycov )
    Pcov = spacetime_attach( p$storage.backend, p$ptr$Pcov )
  }
  if ( exists("TIME", p$variables) ) {
    Ytime = spacetime_attach( p$storage.backend, p$ptr$Ytime )
    Ptime = spacetime_attach( p$storage.backend, p$ptr$Ptime )
  }  

  if ( p$storage.backend != "bigmemory.ram" ) {
    # force copy into RAM to reduce thrashing
    Sloc = Sloc[]
    Yloc = Yloc[]
    Y = Y[]
  }

  if (p$spacetime_engine=="habitat") {
    Ylogit = spacetime_attach( p$storage.backend, p$ptr$Ylogit )
    Plogit = spacetime_attach( p$storage.backend, p$ptr$Plogit )
    Plogitsd = spacetime_attach( p$storage.backend, p$ptr$Plogitsd )
  }

  Yi = spacetime_attach( p$storage.backend, p$ptr$Yi )
  Yi = as.vector(Yi[])  #force copy to RAM as a vector

  # main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    Si = p$runs[ iip, "locs" ]
    if ( is.infinite( Sflag[Si] ) ) next() 
    if ( !is.nan( Sflag[Si] ) ) next() 
    Sflag[Si] = Inf   # over-written below if successful else if a run fails it does not get revisited 
    print( iip )

    # find data nearest S[Si,] and with sufficient data
    dlon = abs( Sloc[Si,1] - Yloc[Yi,1] ) 
    dlat = abs( Sloc[Si,2] - Yloc[Yi,2] ) 
    U =  which( dlon  <= p$dist.max  & dlat <= p$dist.max )
    dist.cur = p$dist.max
    ndata = length(U)
 

    if (ndata > p$n.min ) {
      if (ndata > p$n.max ) {
        Uj = U[ .Internal( sample( ndata, p$n.max, replace=FALSE, prob=NULL)) ]  
      } else {
        Uj = U
      }
      #fast model with gstat to get a crude est of distance scale
      Yyy = data.frame( cbind( p$spacetime_family$linkfun(Y[Uj]), Yloc[Uj,] ))
      names( Yyy) = c("y", "plon", "plat" )
      vEm = try( variogram( y ~ 1, locations=~plon+plat, data=Yyy, cutoff=p$dist.max, width=p$dist.max/15 ) ) # empirical variogram
      if (! inherits(vEm, "try-error")) {
        Yvar = var( Yyy$y, na.rm=TRUE )
        vMod0 = try( vgm(psill=Yvar*0.8, model="Gau", range=p$dist.max, nugget=Yvar*0.2 ) )
        if ( !inherits(vMod0, "try-error")) {
          vFitgs =  try( fit.variogram( vEm, vMod0, fit.sills=TRUE, fit.ranges=TRUE ) ) 
          # plot(vEm, model=vFitgs, add=T)
          if  (! inherits(vFitgs, "try-error") ) {
            dist.cur = max(1, geoR::practicalRange("gaus", phi=vFitgs$range[2] ) )
            U = which( dlon  <= dist.cur  & dlat <= dist.cur )
            ndata =length(U)
          }
        }
        vFit = vMod0 = Yvar = vEm = Yyy = Uj = NULL
      }
    }

    # as a backup .. find data withing a given distance / number 
    if (ndata < p$n.min | ndata > p$n.max | dist.cur < p$dist.min | dist.cur > p$dist.max ) { 
    
      if ( ndata < p$n.min )  {
        for ( usamp in p$upsampling )  {
          dist.cur = p$dist.median * usamp
          U = which( dlon < dist.cur & dlat < dist.cur ) # faster to take a block 
          ndata = length(U)
          if ( ndata >= p$n.min ) {
            if (ndata >= p$n.max) {
              U = U[ .Internal( sample( length(U), p$n.max, replace=FALSE, prob=NULL)) ] 
              ndata = p$n.max
              break()
            }
          }
        }
      } else {
        if ( ndata <= p$n.max * 1.5 ) { # if close to p$n.max, subsample quickly 
          if ( ndata > p$n.max) { 
            U = U[ .Internal( sample( length(U), p$n.max, replace=FALSE, prob=NULL)) ] 
            ndata = p$n.max
          } 
        } else {
          for ( dsamp in p$downsampling )  { # lots of data .. downsample
            dist.cur = p$dist.median * dsamp
            U = which( dlon < dist.cur & dlat < dist.cur )# faster to take a block 
            ndata = length(U)
            if ( ndata <= p$n.max ) break()
            if ( dist.cur <= p$dist.min ) {
              # reached lower limit in distance, taking a subsample instead
              U = which( dlon < p$dist.min & dlat < p$dist.min ) # faster to take a block 
              U = U[ .Internal( sample( length(U), p$n.max, replace=FALSE, prob=NULL)) ]
              ndata = length(U)
              break()
            }
          }
        }
      } 

    }

    rm (dlon, dlat); gc()

    # final check
    ndata = length(U)
    if ((ndata < p$n.min) | (ndata > p$n.max) ) next()
    YiU = Yi[U]  
    # So, YiU and dist_prediction determine the data entering into local model construction
    # dist_model = dist.cur

    dist_prediction = min( p$spacetime_distance_prediction, dist.cur ) # do not predict greater than p$spacetime_distance_prediction

    # construct prediction/output grid area ('pa')
    windowsize.half = floor(dist_prediction/p$pres) # convert distance to discretized increments of row/col indices

    pa_w = -windowsize.half : windowsize.half
    pa_w_n = length(pa_w)
    iwplon = p$rcS[Si,1] + pa_w
    iwplat = p$rcS[Si,2] + pa_w
    pa = NULL
    pa = data.frame( iplon = rep.int(iwplon, pa_w_n) , 
                     iplat = rep.int(iwplat, rep.int(pa_w_n, pa_w_n)) )
    rm(iwplon, iwplat, pa_w)

    bad = which( (pa$iplon < 1 & pa$iplon > p$nplons) | (pa$iplat < 1 & pa$iplat > p$nplats) )
    if (length(bad) > 0 ) pa = pa[-bad,]
    if (nrow(pa)< 5) next()
    
    rc_local = paste(pa$iplon, pa$iplat, sep = "~")
    pa$i = match(rc_local, p$rcP$rc)
    
    bad = which( !is.finite(pa$i))
    if (length(bad) > 0 ) pa = pa[-bad,]

    pa_n = nrow(pa)
    if ( pa_n < 5) next()

      if (0) {
        Sloc = spacetime_attach( p$storage.backend, p$ptr$Sloc )
        Yloc = spacetime_attach( p$storage.backend, p$ptr$Yloc )
        plot( Yloc[U,1]~ Yloc[U,2], col="red", pch=".") # all data
        points( Yloc[YiU,1] ~ Yloc[YiU,2], col="green" )  # with covars and no other data issues
        points( Sloc[Si,1] ~ Sloc[Si,2], col="blue" ) # statistical locations
        points( p$plons[p$rcS[Si,1]] ~ p$plats[p$rcS[Si,2]] , col="purple", pch=25, cex=2 ) # check on p$rcS indexing
        points( p$plons[pa$iplon] ~ p$plats[ pa$iplat] , col="cyan", pch=".", cex=0.01 ) # check on Proc iplat indexing
        points( Ploc[pa$i,1] ~ Ploc[ pa$i, 2] , col="black", pch=20, cex=0.7 ) # check on pa$i indexing -- prediction locations
      }
    rm(rc_local)
   
    pa$plon = Ploc[ pa$i, 1]
    pa$plat = Ploc[ pa$i, 2]

    # prediction covariates i.e., independent variables/ covariates
    pvars = c("plon", "plat", "i")
    if (exists("COV", p$variables)) {
      pvars = c( pvars, p$variables$COV)
      for (ci in 1:length(p$variables$COV)) {
        pa[,p$variables$COV[ci]] = Pcov[ pa$i, ci ]
      }
    }
    pa = pa[, pvars]

    if ( exists("TIME", p$variables) ) {
      pa = cbind( pa[ rep.int(1:pa_n, length(Ptime)), ], 
                       rep.int(Ptime[], rep(pa_n,length(Ptime) )) )
      names(pa) = c( pvars, p$variables$TIME )
      if ( p$variables$TIME != "yr" ) pa$yr = trunc( pa[,p$variables$TIME] )
      if (exists("nw", p)) {
        # where time exists and there are seasonal components, 
        pa$dyear = pa[, p$variables$TIME] - pa$yr  # fractional year
        # additional variables are needed: cos.w, sin.w, etc.. 
        # to add an offset to a trig function (b) must add cos to a sin function
        # y ~ a + c*sin(x+b)
        # y ~ a + c*sin(b)*cos(x) + c*cos(b)*sin(x)  
        #   .. as C*sin(x+b) = C*( cos(b) * sin(x) + sin(b) * cos(x) )
        # y ~ b0 + b1*x1 + b2*x2
        # where: 
        #   a = b0
        #   c^2 = b1^2 + b2^2 = c^2*(sin^2(b) + cos^2(b))
        #   c = sqrt(b1^2 + b2^2)
        #   b1/b2 = tan(b)  
        #   b = arctan(b1/b2)
        pa$cos.w  = cos( pa[,p$variables$TIME] )
        pa$sin.w  = sin( pa[,p$variables$TIME] )
        # compute aditional harmonics only if required (to try to speed things up a bit)
        if ( p$spacetime_engine %in% c( "harmonics.2", "harmonics.3"  ) ) {
          pa$cos.w2 = cos( 2*pa[,p$variables$TIME] )
          pa$sin.w2 = sin( 2*pa[,p$variables$TIME] )
        }
        if ( p$spacetime_engine %in% c( "harmonics.3"  ) ) {
          pa$cos.w3 = cos( 3*pa[,p$variables$TIME] )
          pa$sin.w3 = sin( 3*pa[,p$variables$TIME] )
        }
      }
    }

    # prep dependent data 
    # reconstruct data for modelling (x) and data for prediction purposes (pa)
    x = data.frame( Y[YiU] )
    names(x) = p$variables$Y
    if ( exists("spacetime_family", p) ) {
      x[, p$variables$Y] = p$spacetime_family$linkfun ( x[, p$variables$Y] ) 
    }

    x$plon = Yloc[YiU,1]
    x$plat = Yloc[YiU,2]
    x$Y_wgt = 1 / (( Sloc[Si,1] - x$plat)**2 + (Sloc[Si,2] - x$plon)**2 )# weight data in space: inverse distance squared
    x$Y_wgt[ which( x$Y_wgt < 1e-3 ) ] = 1e-3
    x$Y_wgt[ which( x$Y_wgt > 1 ) ] = 1
    
    if (exists("COV", p$variables)) {
      for (i in 1:length(p$variables$COV )) x[, p$variables$COV[i] ] = Ycov[YiU,i]
    }
     
    if (exists("TIME", p$variables)) {
      x[, p$variables$TIME ] = Ytime[YiU,] 
      if ( p$variables$TIME != "yr" ) x$yr = trunc( x[, p$variables$TIME])
      if (exists("nw", p)) {
        x$dyear = x[, p$variables$TIME] - x$yr
        x$cos.w  = cos( 2*pi*x[,p$variables$TIME] )
        x$sin.w  = sin( 2*pi*x[,p$variables$TIME] )
        if ( p$spacetime_engine %in% c( "harmonics.2", "harmonics.3"  ) ) {
          x$cos.w2 = cos( 2*x[,p$variables$TIME] )
          x$sin.w2 = sin( 2*x[,p$variables$TIME] )
        }
        if ( p$spacetime_engine %in% c( "harmonics.3"  ) ) {
          x$cos.w3 = cos( 3*x[,p$variables$TIME] )
          x$sin.w3 = sin( 3*x[,p$variables$TIME] )
        }
      }
    }

    # model and prediction
    res =NULL
  
    if (p$spacetime_engine=="gam" )  res = spacetime__gam( p, x, pa )
    if (p$spacetime_engine=="habitat") res = spacetime__habitat( p, x, pa )
    if (p$spacetime_engine=="kernel.density") res = spacetime__kerneldensity( p, x, pa )
    if (p$spacetime_engine=="LaplacesDemon") res = spacetime__LaplacesDemon( p, x, pa )
    if (p$spacetime_engine=="inla") res = spacetime__inla( p, x, pa )


    if (0) {
      lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
       
      lattice::levelplot( mean ~ plon + plat, data=res$predictions, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
   
      for( i in sort(unique(res$predictions[,p$variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
    }

    rm(x); gc()
    if ( is.null(res)) next()
   
    if (exists("spacetime_family", p)) {
      res$predictions$mean = p$spacetime_family$linkinv( res$predictions$mean )
      res$predictions$sd  =  p$spacetime_family$linkinv( res$predictions$sd )
    }

    if (exists( "quantile_bounds", p)) {
      tq = quantile( Y[YiU], probs=p$quantile_bounds, na.rm=TRUE  )
      toolow  = which( res$predictions$mean < tq[1] )
      toohigh = which( res$predictions$mean > tq[2] )
      if (length( toolow) > 0)  res$predictions$mean[ toolow] = tq[1]
      if (length( toohigh) > 0) res$predictions$mean[ toohigh] = tq[2]
    }
    
    ii = which( is.finite(res$predictions$mean+res$predictions$sd))
    if (length(ii) < 5) next()  # looks to be a faulty solution


    # stats collator
    if (!exists("spacetime_stats",  res) ) res$spacetime_stats = list()
    if (exists("spacetime_variogram_engine", p) ) {
      sp.stat = try( spacetime_variogram(  xy=Yloc[U,], z=Y[U], methods=p$spacetime_variogram_engine ) )
      if (!is.null(sp.stat) && !inherits(sp.stat, "try-error") ) ){
        res$spacetime_stats["sdSpatial"] = sqrt( sp.stat[[p$spacetime_variogram_engine]]$varSpatial )
        res$spacetime_stats["sdObs"] = sqrt( sp.stat[[p$spacetime_variogram_engine]]$varObs )
        res$spacetime_stats["range"] = sp.stat[[p$spacetime_variogram_engine]]$range
        res$spacetime_stats["phi"] = sp.stat[[p$spacetime_variogram_engine]]$phi
        res$spacetime_stats["nu"] = sp.stat[[p$spacetime_variogram_engine]]$nu
      } 
      sp.stat = NULL
    }
  
    if ( exists("TIME", p$variables) ){
      # annual ts, seasonally centered and spatially 
      # pa_i = which( Sloc[Si,1]==Ploc[,1] & Sloc[Si,2]==Ploc[,2] )
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
          if (!is.null(ts.stat) && !inherits(ts.stat, "try-error") ) {
            res$spacetime_stats["ar_timerange"] = ts.stat$quantilePeriod 
            if (length(which (is.finite(pac$mean))) > 5 ) {
              ar1 = try( lm( pac$mean[1:(length(piid) - 1)] ~ pac$mean[2:(length(piid))] + 0, na.action="na.omit") )
              if (!inherits(ts.stat, "try-error")) res$spacetime_stats["ar_1"] = coef( ar1 )
        } } } 
        rm ( pac, piid )
      } 
      rm(pac_i)
    }
       

    # update SD estimates of predictions with those from other locations via the
    # incremental  method ("online algorithm") of mean estimation after Knuth ;
    # see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    # update means: inverse-variance weighting   
    # see https://en.wikipedia.org/wiki/Inverse-variance_weighting
   
    npred = nrow(res$predictions)

    if ( ! exists("TIME", p$variables) ) {

      u = which( is.finite( P[res$predictions$i] ) )  # these have data already .. update
      if ( length( u ) > 1 ) {
        ui = res$predictions$i[u]  # locations of P to modify
        Pn[ui] = Pn[ui] + 1 # update counts
        stdev_update =  Psd[ui] + ( res$predictions$sd[u] -  Psd[ui] ) / Pn[ui]
        means_update = ( P[ui] / Psd[ui]^2 + res$predictions$mean[u] / res$predictions$sd[u]^2 ) / ( Psd[ui]^(-2) + res$predictions$sd[u]^(-2) )
        mm = which(is.finite( means_update + stdev_update ))
        if( length(mm)> 0) {
          iumm = ui[mm]
          Psd[iumm] = stdev_update[mm]
          P  [iumm] = means_update[mm]
        } 
        stdev_update = NULL
        means_update = NULL

        if (p$spacetime_engine=="habitat") {
          logit_stdev_update =  Plogitsd[ui] + ( res$predictions$logitsd[u] -  Plogitsd[ui] ) / Pn[ui]
          logit_means_update = ( Plogit[ui] / Plogitsd[ui]^2 + res$predictions$logitmean[u] / res$predictions$logitsd[u]^2 ) / ( Plogitsd[ui]^(-2) + res$predictions$logitsd[u]^(-2) )
          mm = which(is.finite( logit_means_update + logit_stdev_update ))
          if( length(mm)> 0) {
            iumm = ui[mm]
            Plogitsd[iumm] = logit_stdev_update[mm]
            Plogit  [iumm] = logit_means_update[mm]
          }
          logit_stdev_update = NULL
          logit_means_update = NULL
        }
        rm(ui, mm, iumm)
      }

      # first time # no data yet
      v = setdiff(1:npred, u)         
      if ( length(v) > 0 ) {
        vi = res$predictions$i[v]
        Pn [vi] = 1
        P  [vi] = res$predictions$mean[v]
        Psd[vi] = res$predictions$sd[v]
        if (p$spacetime_engine=="habitat") {
          Plogit  [vi] = res$predictions$logitmean[v]
          Plogitsd[vi] = res$predictions$logitsd[v]
        }
      }
    }

    if ( exists("TIME", p$variables) ) {
      u = which( is.finite( P[res$predictions$i,1] ) )  # these have data already .. update
      nu = length( u ) 
      if ( nu > 1 ) {  # ignore if only one point .. mostly because it can cause issues with matrix form .. 
        # locations of P to modify
        ui = sort(unique(res$predictions$i[u]))
        nc = ncol(P)
        if (p$storage.backend == "ff" ) {
          add.ff(Pn, 1, ui, 1:nc ) # same as Pn[ui,] = Pn[ui]+1 but 2X faster
        } else {
          Pn[ui,] = Pn[ui,] + 1
        }
        stdev_update =  Psd[ui,] + ( res$predictions$sd[u] -  Psd[ui,] ) / Pn[ui,]
        means_update = ( P[ui,] / Psd[ui,]^2 + res$predictions$mean[u] / res$predictions$sd[u]^2 ) / 
          ( Psd[ui,]^(-2) + res$predictions$sd[u]^(-2) )
        
        updates = means_update + stdev_update 
        if (!is.matrix(updates)) next()

        mm = which( is.finite( rowSums(updates)))  # created when preds go outside quantile bounds .. this removes all data from a given location rather than the space-time .. severe but likely due to a poor prediction and so remove all (it is also faster this way as few manipulations)
        if( length(mm)> 0) {
          iumm = ui[mm] 
          Psd[iumm,] = stdev_update[mm,]
          P  [iumm,] = means_update[mm,]
          iumm = NULL
        } 
        stdev_update = NULL
        means_update = NULL
        if (p$spacetime_engine=="habitat") {
          logit_stdev_update =  Plogitsd[ui,] + ( res$predictions$logitsd[u] -  Plogitsd[ui,] ) / Pn[ui]
          logit_means_update = ( Plogit[ui,] / Plogitsd[ui,]^2 + res$predictions$logitmean[u] / res$predictions$logitsd[u]^2 ) / ( Plogitsd[ui,]^(-2) + res$predictions$logitsd[u]^(-2) )
          updates = logit_means_update + logit_stdev_update
          if (!is.matrix(updates)) next()
          mm = which( is.finite( rowSums(updates)))  # created when preds go outside quantile bounds .. this removes 
          if( length(mm)> 0) {
            iumm = ui[mm]
            Plogitsd[iumm,] = logit_stdev_update[mm,]
            Plogit  [iumm,] = logit_means_update[mm,]
            iumm = NULL
          } 
          logit_stdev_update = NULL
          logit_means_update = NULL
        }
        rm(ui, mm)

      }

      # do this as a second pass in case NA's were introduced by the update .. unlikely , but just in case
      v = which( !is.finite( P[res$predictions$i,1] ) )  # these have data already .. update
      nv = length(v)          # no data yet
      if ( nv > 0 ) {
        vi = sort(unique(res$predictions$i[v]))
        Pn [vi,] = 1
        P  [vi,] = res$predictions$mean[v]
        Psd[vi,] = res$predictions$sd[v]
        if (p$spacetime_engine=="habitat") {
          Plogit  [vi,] = res$predictions$logitmean[v]
          Plogitsd[vi,] = res$predictions$logitsd[v]
        }
        rm(vi)
      } 
    }


    # save stats
    for ( k in 1: length(p$statsvars) ) {
      if (exists( p$statsvars[k], res$spacetime_stats )) {
        S[Si,k] = res$spacetime_stats[[ p$statsvars[k] ]]
      }
    }
    
      if (0) {
     
        lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
        
        for( i in sort(unique(res$predictions[,p$variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
      
        for (i in 1:p$nt) {
          print( lattice::levelplot( P[pa$i,i] ~ Ploc[pa$i,1] + Ploc[ pa$i, 2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
        }

        v = res$predictions
        if ( exists("TIME", p$variables) ){
          v = v[which( v[,p$variables$TIME]==1990.55),]
        }
        require(lattice)
        levelplot( mean ~ plon+plat, v, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
      }
   
    res = NULL
    pa = NULL

    # ----------------------
    # do last. it is an indicator of completion of all tares$predictionssks .. restarts would be broken otherwise
    Sflag[Si] = 1  # done .. any finite value would do

  }  # end for loop
  
  invisible()

}

