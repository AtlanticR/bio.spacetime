
spacetime_interpolate = function( ip=NULL, p ) {

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  #---------------------
  # data for modelling
  # dependent vars # already link-transformed in spacetime_db("dependent")
  Y = attach.big.matrix( p$ptr$Y )
  Yloc = attach.big.matrix( p$ptr$Yloc )

  P = attach.big.matrix( p$ptr$P )
  Pn = attach.big.matrix( p$ptr$Pn )
  Psd = attach.big.matrix( p$ptr$Psd )
  Ploc = attach.big.matrix( p$ptr$Ploc )
  
  Sloc = attach.big.matrix( p$ptr$Sloc ) # statistical output locations
  S = attach.big.matrix( p$ptr$S )  # statistical outputs inside loop to safely save data and pass onto other processes

  Yi = 1:length(Y)
  bad = which( !is.finite( Y[]))
  if (length(bad)> 0 ) Yi[bad] = NA

  # data locations
  bad = which( !is.finite( rowSums(Yloc[])))
  if (length(bad)> 0 ) Yi[bad] = NA

# data locations
  if (exists("COV", p$variables)) {
    Ycov = attach.big.matrix( p$ptr$Ycov )
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
    Ytime = attach.big.matrix( p$ptr$Ytime )
    bad = which( !is.finite( Ytime[] ))
    if (length(bad)> 0 ) Yi[bad] = NA
    Yi = na.omit(Yi)
  }

  #---------------------
  # prediction locations and covariates

  rcP = data.frame( cbind( 
    Prow = (Ploc[,1]-p$plons[1])/p$pres + 1,  
    Pcol = (Ploc[,2]-p$plats[1])/p$pres + 1) )
  # rcP$i =1:nrow(rcP) # row index
  rcP$rc = paste( rcP$Prow, rcP$Pcol, sep="~")
  rcP$Prow = rcP$Pcol = NULL
  gc()

  #-----------------
  # row, col indices
  
  rcS = data.frame( cbind( 
    Srow = (Sloc[,1]-p$plons[1])/p$pres + 1,  
    Scol = (Sloc[,2]-p$plats[1])/p$pres + 1))


  # main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    Si = p$runs[ iip, "locs" ]
    if ( is.infinite( S[Si,1] ) ) next() 
    if ( !is.nan( S[Si,1] ) ) next() 
    S[Si,1] = Inf   # over-written below if successful else if a run fails it does not get revisited 
 
    # find data withing a given distance / number 
    pib = point_in_block( Sloc=Sloc, Si=Si, Yloc=Yloc, Yi=Yi, 
      dist.max=p$dist.max, dist.min=p$dist.min, 
      n.min=p$n.min, n.max=p$n.max, 
      upsampling=p$upsampling, downsampling=p$downsampling, resize=TRUE ) 
    if ( is.null(pib)) {
      next()
    } else {
      dist.cur = pib$dist
      U = pib$U
      rm(pib); gc()
    }
    ndata = length(U)
    if ((ndata < p$n.min) | (ndata > p$n.max) ) next()
    YiU = Yi[U]  

    # So, YiU and dist.cur determine the data entering into local model construction
    # but predictions are expensive, esp in space-time and so a smaller spatial raange is selected
    # this is determined by p$spacetime.prediction.dist.min
    dist_prediction = min( dist.cur*p$spacetime.prediction.range.proportion,  p$spacetime.prediction.dist.min )
    windowsize.half = floor(dist_prediction/p$pres) # convert distance to discretized increments of row/col indices
    pa_offsets = -windowsize.half : windowsize.half
    pa_offsets_n = length(pa_offsets)
    pa_prow = rcS[Si,1] + pa_offsets
    pa_pcol = rcS[Si,2] + pa_offsets

    pa = data.frame( Prow = rep.int(pa_prow, pa_offsets_n) , 
                     Pcol = rep.int(pa_pcol, rep.int(pa_offsets_n,pa_offsets_n)) )

    bad = which( (pa$Prow < 1 & pa$Prow > p$nplons) | (pa$Pcol < 1 & pa$Pcol > p$nplats) )
    if (length(bad) > 0 ) pa = pa[-bad,]
    if (nrow(pa)< 5) next()

    # 
    pc_rc = paste( pa$Prow, pa$Pcol, sep="~" )
    pa$i = match( pc_rc, rcP$rc)
    bad = which( !is.finite(pa$i))
    if (length(bad) > 0 ) pa = pa[-bad,]
    if (nrow(pa)< 5) next()
    pa$plon = Ploc[ pa$i, 1]
    pa$plat = Ploc[ pa$i, 2]

    # prediction covariates i.e., independent variables/ covariates
    if (exists("COV", p$variables)) {
      Pcov = attach.big.matrix( p$ptr$Pcov )
      for (ci in 1:length(p$variables$COV)) {
        pa[,p$variables$COV[ci]] = Pcov[ pa$i, ci ]
      }
    }

    if (0) {
      plot( Yloc[U,1]~ Yloc[U,2], col="red", pch=".") # all data
      points( Yloc[YiU,1] ~ Yloc[YiU,2], col="green" )  # with covars and no other data issues
      points( Sloc[Si,1] ~ Sloc[Si,2], col="blue" ) # statistical locations
      points( p$plons[rcS[Si,1]] ~ p$plats[rcS[Si,2]] , col="purple", pch=25, cex=2 ) # check on rcS indexing
      points( p$plons[pa$Prow] ~ p$plats[ pa$Pcol] , col="cyan", pch=".", cex=0.01 ) # check on Proc Pcol indexing
      points( Ploc[pa$i,1] ~ Ploc[ pa$i, 2] , col="black", pch=20, cex=0.7 ) # check on pa$i indexing -- prediction locations
    }

  # ----------------------
  # prediction
    Pi = 1:nrow( Ploc ) 
    pa$i = Pi[ pa$i ]  ## flag to ignore
    kP = na.omit(pa$i)
    if ( length( kP) < 5 ) next()
    pa = pa[ which(is.finite( pa$i)), ] 
    pa_n = nrow(pa)

    # construct prediction/output grid 
    if ( ! exists("TIME", p$variables) ) {
      if (  exists("COV", p$variables) ) {
        newdata = pa[ , c("plon", "plat", "i", p$variables$COV)]
      } else {
        newdata = pa[ , c("plon", "plat", "i")]
      }
   } else {
      Ptime = attach.big.matrix( p$ptr$Ptime )
      if ( exists("COV", p$variables) ) {
        pa_t = pa[ , c("plon", "plat", "i", p$variables$COV)]
        newdata = cbind( pa_t[ rep.int(1:pa_n, length(Ptime)), ], rep.int(Ptime[], pa_n ))
        names(newdata) = c("plon", "plat", "i", p$variables$COV, "tiyr" )
      } else {
        pa_t = pa[ , c("plon", "plat", "i")]
        newdata = cbind( pa_t[ rep.int(1:pa_n, length(Ptime)), ], rep.int(Ptime[], pa_n ))
        names(newdata) = c("plon", "plat", "i", "tiyr" )
      }
    }


    if ( p$spacetime_engine %in% 
      c( "harmonics.1", "harmonics.2", "harmonics.3", "harmonics.1.depth",
         "seasonal.basic", "seasonal.smoothed", "annual", "gam"  ) ) {
      res = spacetime__harmonics( p, YiU, Si, newdata )
      if ( is.null(res)) next()     

      if (exists("variogram.engine", p) ) {
        sp.stat = spacetime_variogram(  Yloc[YiU,], Y[YiU], methods=p$spacetime_variogram_engine )
        res$spacetime_stats["sdSpatial"] = sqrt( sp.stat[[p$variogram.engine]]$varSpatial )
        res$spacetime_stats["sdObs"] = sqrt( sp.stat[[p$variogram.engine]]$varObs )
        res$spacetime_stats["range"] = sp.stat[[p$variogram.engine]]$range
        res$spacetime_stats["phi"] = sp.stat[[p$variogram.engine]]$phi
        res$spacetime_stats["nu"] = sp.stat[[p$variogram.engine]]$nu
      }

      if ( exists("TIME", p$variables) ){
        if (exists("timeseries.engine", p) ) {
           # annual ts, seasonally centered and spatially 
          pa_S = paste( rcS[Si,1], rcS[Si,2], sep="~" )
          pa_i = match( pa_S, rcP$rc)
          if (length(pa_i)==1) {

            pii = which(pa$i== pa_i)
            pac = pa[pii,]
            pac$dyr = pac[, p$variables$TIME] - trunc(pac[, p$variables$TIME] )
            piid = which( zapsmall( pac$dyr - p$dyear_centre) == 00 )
            pac = pac[ piid, c(p$variables$TIME, "mean")]
            pac = pac[ order(pac[,p$variables$TIME]),]
            if (length(pii) > 5 ) {
              ts.stat = spacetime_timeseries( pac$mean, method="fft" )
              res$spacetime_stats["ar_timerange"] = ts.stat$quantilePeriod 
              res$spacetime_stats["ar_1"] = coef( lm( pac$mean[1:(length(piid) - 1)] ~ pac$mean[2:(length(piid))] + 
    0 ) )
       
            }
          }
        }
      }
    }

    if (p$spacetime_engine=="kernel.density") res = spacetime__kerneldensity(  )
    if (p$spacetime_engine=="LaplacesDemon") res = spacetime__LaplacesDemon(  )
    if (p$spacetime_engine=="inla") res = spacetime__inla()

    if ( is.null(res)) next()

    pa = merge( res$predictions[,c("i", "mean", "sd")], pa[,c("i", "Prow", "Pcol")], by="i", all.x=TRUE, all.y=FALSE, sort=FALSE )
    
    spacetime_stats = res$spacetime_stats
    
    rm(res, newdata); gc()

    if (debugrun) {
      pa = merge( res$predictions[,c("i", "mean", "sd", "tiyr", "plon","plat")], pa[,c("i", "Prow", "Pcol")], by="i", all.x=TRUE, all.y=FALSE, sort=FALSE )
      it = which(pa$tiyr==1990.55)
      x11(); levelplot( mean ~ plon+plat, pa[it,], aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
      it = which(pa$tiyr==2010.55)
      x11(); levelplot( mean   ~ plon+plat, pa[it,], aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }

    good = which( is.finite( rowSums(pa) ) )
    if (length(good) < 1) next()
    pa = pa[good,]

    # ----------------------
    # update P (predictions)
    # maybe faster to make copy first and then do a rowSums ?
    test = rowSums( P[pa$i,] ) 
    u = which( is.finite( test ) )  # these have data already .. update
    if ( length( u ) > 0 ) {
      ui = pa$i[u]  # locations of P to modify
      Pn[ui,] = Pn[ui,] + 1 # update counts
      # update SD estimates of predictions with those from other locations via the
      # incremental  method ("online algorithm") of mean estimation after Knuth ;
      # see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
      stdev_update =  Psd[ui,] + ( pa$sd[u] -  Psd[ui,] ) / Pn[ui,]
      # update means: inverse-variance weighting   https://en.wikipedia.org/wiki/Inverse-variance_weighting
      means_update = ( P[ui,] / Psd[ui,]^2 + pa$mean[u] / pa$sd[u]^2 ) / ( Psd[ui,]^(-2) + pa$sd[u]^(-2) )
      mm = which(is.finite( means_update + stdev_update ))
      if( length(mm)> 0) {
        # actual updates occur after everything has been computed first
        iumm = ui[mm]
        Psd[iumm,] = stdev_update[mm]
        P  [iumm,] = means_update[mm]
      }
    }

    # do this as a second pass in case NA's were introduced by the update .. unlikely , but just in case
    test = rowSums( P[pa$i,] )
    f = which( !is.finite( test ) ) # first time
    if ( length(f) > 0 ) {
      fi = pa$i[f]
      Pn [fi,] = 1
      P  [fi,] = pa$mean[f]
      Psd[fi,] = pa$sd[f]
    }
  
    rm( pa ) ; gc()

    #########

    # print( "Saving summary statisitics" )
    # save statistics last as this is an indicator of completion of all tasks .. restarts would be broken otherwise
    for ( k in 1: length(p$statsvars) ) {
      if (exists( p$statsvars[k], spacetime_stats )) S[Si,k] = spacetime_stats[[ p$statsvars[k] ]]
    }

  }  # end for loop

}



