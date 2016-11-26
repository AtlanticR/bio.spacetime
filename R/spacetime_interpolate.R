
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

  # misc intermediate calcs to be done outside of parallel loops
  upsampling = sort( p$sampling[ which( p$sampling > 1 ) ] )
  upsampling = upsampling[ which(upsampling*p$spacetime_distance_scale <= p$spacetime_distance_max )]
  downsampling = sort( p$sampling[ which( p$sampling < 1) ] , decreasing=TRUE )
  downsampling = downsampling[ which(downsampling*p$spacetime_distance_scale >= p$spacetime_distance_min )]

  # for 2D methods, treat time as independent timeslices
  if ( exists("TIME", p$variables)) {
    p$ts = Ptime[]
  } else {
    p$ts = 1
  }

  # used by "fields":
  theta.grid = 10^seq( -6, 6, by=0.5) * p$spacetime_distance_scale # maxdist is aprox magnitude of the phi parameter
  lambda.grid = 10^seq( -9, 3, by=0.5) 

  #-----------------
  # row, col indices
  # statistical output locations
  rcS = data.frame( cbind( 
    Srow = (Sloc[,1]-p$plons[1])/p$pres + 1,  
    Scol = (Sloc[,2]-p$plats[1])/p$pres + 1))

  #---------------------
  # prediction locations and covariates
  rcP = data.frame( cbind( 
    Prow = (Ploc[,1]-p$plons[1])/p$pres + 1,  
    Pcol = (Ploc[,2]-p$plats[1])/p$pres + 1) )
  rcPid = paste( rcP$Prow, rcP$Pcol, sep="~")
  rm(rcP)

  gc()


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
    U =  which( dlon  <= p$spacetime_distance_scale  & dlat <= p$spacetime_distance_scale )
    spacetime_distance_cur = p$spacetime_distance_scale
    ndata = length(U)
  
    o = ores = NULL

    if (ndata > p$n.min ) {
      if (ndata > p$n.max ) {
        Uj = U[ .Internal( sample( ndata, p$n.max, replace=FALSE, prob=NULL)) ]  
      } else {
        Uj = U
      }
      o = try( spacetime_variogram( xy=Yloc[Uj,], z=p$spacetime_family$linkfun(Y[Uj]), methods=p$spacetime_engine.variogram ) )
      if (!inherits(o, "try-error")) {
        if ( !is.null(o)) {
          spacetime_distance_cur = min( max(1, o[[p$spacetime_engine.variogram]][["range"]] ), p$spacetime_distance_scale ) 
          U = which( dlon  <= spacetime_distance_cur  & dlat <= spacetime_distance_cur )
          ndata =length(U)
          smoothness0 = o[[p$spacetime_engine.variogram]][["nu"]]
          ores = o[[p$spacetime_engine.variogram]] # store current best estimate of variogram characteristics
          if (0) {
            if (p$spacetime_engine.variogram == "fast" ) { 
              plot(o[[p$spacetime_engine.variogram]][["vgm"]], 
                model=RMmatern( nu=o$fast$nu, var=o$fast$varSpatial, scale=o$fast$phi * (sqrt(o$fast$nu*2) )) + RMnugget(var=o$fast$varObs) )
            }
          }
        }   
      }
    }

    # if insufficient data found within the "range" fall back to a brute force search until criteria are met
    if (ndata < p$n.min | ndata > p$n.max | spacetime_distance_cur < p$spacetime_distance_min | spacetime_distance_cur > p$spacetime_distance_max ) { 
      if ( ndata < p$n.min )  {
        for ( usamp in upsampling )  {
          spacetime_distance_cur = p$spacetime_distance_scale * usamp
          U = which( dlon < spacetime_distance_cur & dlat < spacetime_distance_cur ) # faster to take a block 
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
          for ( dsamp in downsampling )  { # lots of data .. downsample
            spacetime_distance_cur = p$spacetime_distance_scale * dsamp
            U = which( dlon < spacetime_distance_cur & dlat < spacetime_distance_cur )# faster to take a block 
            ndata = length(U)
            if ( ndata <= p$n.max ) break()
            if ( spacetime_distance_cur <= p$spacetime_distance_min ) {
              # reached lower limit in distance, taking a subsample instead
              U = which( dlon < p$spacetime_distance_min & dlat < p$spacetime_distance_min ) # faster to take a block 
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
    # dist_model = spacetime_distance_cur

    dist_prediction = min( p$spacetime_distance_prediction, spacetime_distance_cur ) # do not predict greater than p$spacetime_distance_prediction

    # construct prediction/output grid area ('pa')
    windowsize.half = floor(dist_prediction/p$pres) # convert distance to discretized increments of row/col indices

    pa_w = -windowsize.half : windowsize.half
    pa_w_n = length(pa_w)
    iwplon = rcS[Si,1] + pa_w
    iwplat = rcS[Si,2] + pa_w
    pa = NULL
    pa = data.frame( iplon = rep.int(iwplon, pa_w_n) , 
                     iplat = rep.int(iwplat, rep.int(pa_w_n, pa_w_n)) )
    rm(iwplon, iwplat, pa_w)

    bad = which( (pa$iplon < 1 & pa$iplon > p$nplons) | (pa$iplat < 1 & pa$iplat > p$nplats) )
    if (length(bad) > 0 ) pa = pa[-bad,]
    if (nrow(pa)< 5) next()
    
    rc_local = paste(pa$iplon, pa$iplat, sep = "~")
    pa$i = match(rc_local, rcPid)
    
    bad = which( !is.finite(pa$i))
    if (length(bad) > 0 ) pa = pa[-bad,]

    pa_n = nrow(pa)
    if ( pa_n < 5) next()

      if (0) {
        # check that position indices are working properly
        Sloc = spacetime_attach( p$storage.backend, p$ptr$Sloc )
        Yloc = spacetime_attach( p$storage.backend, p$ptr$Yloc )
        plot( Yloc[U,1]~ Yloc[U,2], col="red", pch=".") # all data
        points( Yloc[YiU,1] ~ Yloc[YiU,2], col="green" )  # with covars and no other data issues
        points( Sloc[Si,1] ~ Sloc[Si,2], col="blue" ) # statistical locations
        points( p$plons[rcS[Si,1]] ~ p$plats[rcS[Si,2]] , col="purple", pch=25, cex=2 ) # check on rcS indexing
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
      # where time exists and there are seasonal components, 
      # additional variables are created/needed here: cos.w, sin.w, etc.. 
      # for harmonic analysis: to add an offset to a trig function (b) must add cos to a sin function
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
      if ("dyear" %in% p$variables$ALL)  pa$dyear = pa[, p$variables$TIME] - pa$yr  # fractional year
      if ("cos.w" %in% p$variables$ALL)  pa$cos.w  = cos( pa[,p$variables$TIME] )
      if ("sin.w" %in% p$variables$ALL)  pa$sin.w  = sin( pa[,p$variables$TIME] )
      if ("cos.w2" %in% p$variables$ALL) pa$cos.w2 = cos( 2*pa[,p$variables$TIME] )
      if ("sin.w2" %in% p$variables$ALL) pa$sin.w2 = sin( 2*pa[,p$variables$TIME] )
      if ("cos.w3" %in% p$variables$ALL) pa$cos.w3 = cos( 3*pa[,p$variables$TIME] )
      if ("sin.w3" %in% p$variables$ALL) pa$sin.w3 = sin( 3*pa[,p$variables$TIME] )
      # more than 3 harmonics would not be advisable .. but you would add them here..
    }
    
    # prep dependent data 
    # reconstruct data for modelling (x) and data for prediction purposes (pa)
    x = data.frame( Y[YiU] )
    names(x) = p$variables$Y
    x[, p$variables$Y] = p$spacetime_family$linkfun ( x[, p$variables$Y] ) 
    if (p$spacetime_engine=="habitat") {
      x[, p$variables$Ylogit ] = p$spacetime_family_logit$linkfun ( x[, p$variables$Ylogit] ) ### -- need to conform with data structure ... check once ready
    }
    x$plon = Yloc[YiU,1]
    x$plat = Yloc[YiU,2]
    x$weights = 1 / (( Sloc[Si,1] - x$plat)**2 + (Sloc[Si,2] - x$plon)**2 )# weight data in space: inverse distance squared
    x$weights[ which( x$weights < 1e-3 ) ] = 1e-3
    x$weights[ which( x$weights > 1 ) ] = 1
    
    if (exists("COV", p$variables)) {
      for (i in 1:length(p$variables$COV )) x[, p$variables$COV[i] ] = Ycov[YiU,i]
    }
     
    if (exists("TIME", p$variables)) {
      x[, p$variables$TIME ] = Ytime[YiU,] 
      if ( p$variables$TIME != "yr" ) x$yr = trunc( x[, p$variables$TIME]) 
      if ("dyear" %in% p$variables$ALL)  x$dyear = x[, p$variables$TIME] - x$yr
      if ("cos.w" %in% p$variables$ALL)  x$cos.w  = cos( 2*pi*x[,p$variables$TIME] )
      if ("sin.w" %in% p$variables$ALL)  x$sin.w  = sin( 2*pi*x[,p$variables$TIME] )
      if ("cos.w2" %in% p$variables$ALL) x$cos.w2 = cos( 2*x[,p$variables$TIME] )
      if ("sin.w2" %in% p$variables$ALL) x$sin.w2 = sin( 2*x[,p$variables$TIME] )
      if ("cos.w3" %in% p$variables$ALL) x$cos.w3 = cos( 3*x[,p$variables$TIME] )
      if ("sin.w3" %in% p$variables$ALL) x$sin.w3 = sin( 3*x[,p$variables$TIME] )
    }

    o = NULL
    o = try( spacetime_variogram( xy=Yloc[U,], z=p$spacetime_family$linkfun(Y[U]), methods=p$spacetime_engine.variogram) )
      if (!inherits(o, "try-error")) {
        if ( !is.null(o) ) {
          if ( exists( p$spacetime_engine.variogram, o )) {
            ores = o[[p$spacetime_engine.variogram]]  # replace with this "tweaked variogram estimate"    
          }  
        }
      } 
      
    smoothness = smoothness0
    if (!is.null(ores)) {
      if ( exists("nu", ores) ) {
        smoothness = ores$nu
      } 
    }

    # model and prediction
    # the following permits user-defined models (might want to use compiler::cmpfun )
    gc()
    res =NULL
    res = switch( p$spacetime_engine, 
      gam = spacetime__gam( p, x, pa ), 
      habitat = spacetime__habitat( p, x, pa ), 
      kernel.density = uspacetime__kerneldensity( p, x, pa, smoothness ),
      bayesx = spacetime__bayesx( p, x, pa ),
      LaplacesDemon = spacetime__LaplacesDemon( p, x, pa ),
      gaussianprocess2Dt = spacetime__gaussianprocess2Dt( p, x, pa ), 
      gaussianprocess = spacetime__gaussianprocess( p, x, pa ), 
      inla = spacetime__inla( p, x, pa ),
      spacetime_engine_user_defined = p$spacetime_engine_user_defined( p, x, pa)
    )
    
    if (0) {
      lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
       
      lattice::levelplot( mean ~ plon + plat, data=res$predictions, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
   
      for( i in sort(unique(res$predictions[,p$variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
    }

    rm(x); gc()
    if ( is.null(res)) next()
   
    res$predictions$mean = p$spacetime_family$linkinv( res$predictions$mean )
    res$predictions$sd   = p$spacetime_family$linkinv( res$predictions$sd )
    if (p$spacetime_engine=="habitat") {
      res$predictions$logitmean = p$spacetime_family_logit$linkinv( res$predictions$logitmean )
      res$predictions$logitsd   = p$spacetime_family_logit$linkinv( res$predictions$logitsd )
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
    
    if (!exists("sdSpatial", res$spacetime_stats)) {
      # some methods can generate spatial stats simultaneously .. 
      # it is faster to keep them all together instead of repeating here
      # field and RandomFields gaussian processes seem most promising ... 
      # default to fields for speed:
      res$spacetime_stats["sdSpatial"] = NA 
      res$spacetime_stats["sdObs"] = NA 
      res$spacetime_stats["range"] = NA
      res$spacetime_stats["phi"] = NA
      res$spacetime_stats["nu"] = NA
      if ( !is.null(ores)) {
        res$spacetime_stats["sdSpatial"] = sqrt( ores[["varSpatial"]] ) 
        res$spacetime_stats["sdObs"] = sqrt(ores[["varObs"]]) 
        res$spacetime_stats["range"] = ores[["range"]]
        res$spacetime_stats["phi"] = ores[["phi"]]
        res$spacetime_stats["nu"] = ores[["nu"]]
      } 
    }
    
    if ( exists("TIME", p$variables) ){
      # annual ts, seasonally centered and spatially 
      # pa_i = which( Sloc[Si,1]==Ploc[,1] & Sloc[Si,2]==Ploc[,2] )
      pac_i = which( res$predictions$plon==Sloc[Si,1] & res$predictions$plat==Sloc[Si,2] )
      # plot( mean~tiyr, res$predictions[pac_i,])
      # plot( mean~tiyr, res$predictions, pch="." )
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
              ar1 = NULL
              ar1 = try( ar( pac$mean, order.max=1 ) )
              if (!inherits(ar1, "try-error")) {
                res$spacetime_stats["ar_1"] = ar1$ar 
              } else {
                ar1 = try( cor( pac$mean[1:(length(piid) - 1)], pac$mean[2:(length(piid))], na.rm=TRUE ) )
                if (!inherits(ar1, "try-error")) res$spacetime_stats["ar_1"] = ar1 
              }
            } 
          } 

          ### Do the logistic model here ! -- if not already done ..
          if (!exists("ts_K", res$spacetime_stats)) {
            # model as a logistic with ts_r, ts_K, etc .. as stats outputs
            
          } 

        } 
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
      u_n = length( u ) 
      if ( u_n > 1 ) {  # ignore if only one point .. mostly because it can cause issues with matrix form .. 
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

