

spacetime.interpolate.xyts = function( ip, p ) {

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  #---------------------
  # data for modelling
  # dependent vars # already link-transformed in spacetime.db("dependent")
  Y = p$ff$Y
  Yi = 1:length(Y)
  bad = which( !is.finite( Y[]))
  if (length(bad)> 0 ) Yi[bad] = NA

  # data locations
  Yloc = p$ff$Yloc
  bad = which( !is.finite( rowSums(Yloc[])))
  if (length(bad)> 0 ) Yi[bad] = NA

  # data locations
  Ytime = p$ff$Ytime
  bad = which( !is.finite( rowSums(Ytime[])))
  if (length(bad)> 0 ) Yi[bad] = NA
  Yi = na.omit(Yi)

  #---------------------
  # prediction locations and covariates
  Ploc = p$ff$Ploc  # prediction locations
  Ptime = p$ff$Ptime  # prediction times
  Pcov = p$ff$Pcov
  
  rcP = data.frame( cbind( 
    Prow = (Ploc[,1]-p$plons[1])/p$pres + 1,  
    Pcol = (Ploc[,2]-p$plats[1])/p$pres + 1) )
  # rcP$i =1:nrow(rcP) # row index
  rcP$rc = paste( rcP$Prow, rcP$Pcol, sep="~")
  rcP$Prow = rcP$Pcol = NULL
  gc()

  #-----------------
  # row, col indices
  Sloc = p$ff$Sloc  # statistical output locations
  rcS = data.frame( cbind( 
    Srow = (Sloc[,1]-p$plons[1])/p$pres + 1,  
    Scol = (Sloc[,2]-p$plats[1])/p$pres + 1))

# choose model formula for GAM-based models
  mf = switch( p$tsmethod ,
    annual = ' t ~ s(yr) ',
    seasonal.basic = ' t ~ s(yr) + s(dyear, bs="cc") ', 
    seasonal.smoothed = ' t ~ s(yr, dyear) + s(yr) + s(dyear, bs="cc")  ', 
    harmonics.1 = ' t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w)  ', 
    harmonics.2 = ' t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) ' , 
    harmonics.3 = ' t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) + s(yr, cos.w3) + s(yr, sin.w3)  + s(cos.w3) + s( sin.w3 ) '
  )

  mf = formula(mf)


  # main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    Si = p$runs[ iip, "locs" ]
    S = p$ff$S  # statistical outputs inside loop to safely save data and pass onto other processes
    if ( is.nan( S[Si,1] ) ) next()
    if ( !is.na( S[Si,1] ) ) next()
    S[Si,1] = NaN   # over-written below if successful else if a run fails it does not get revisited 
    close(S) # close right away such that other processes can write to S
    
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

    #  NOTE:: by default, all areas chosen to predict within the window.. but if covariates are involved,
    #  this can be done only where covariates exists .. so next step is to determine this and predict
    #  over the correct area.
    #  Once all predictions are complete, simple (kernal-based?) interpolation
    #  for areas without covariates can be completed
    
    windowsize.half = floor(dist.cur/p$pres) # convert distance to discretized increments of row/col indices
    pa_offsets = -windowsize.half : windowsize.half
    pa = expand.grid( Prow = rcS[Si,1] + pa_offsets, Pcol = rcS[Si,2] + pa_offsets ) # row,col coords
    bad = which( (pa$Prow < 1 & pa$Prow > p$nplons) | (pa$Pcol < 1 & pa$Pcol > p$nplats) )
    if (length(bad) > 0 ) pa = pa[-bad,]
    if (nrow(pa)< p$n.min) next()
    pc_rc = paste( pa$Prow, pa$Pcol, sep="~" )
    pa$i = match( pc_rc, rcP$rc)
    bad = which( !is.finite(pa$i))
    if (length(bad) > 0 ) pa = pa[-bad,]
    if (nrow(pa)< p$n.min) next()
    pa$plon = Ploc[ pa$i, 1]
    pa$plat = Ploc[ pa$i, 2]
 
    if (0) {
      plot( Yloc[U,1]~ Yloc[U,2], col="red", pch=".") # all data
      points( Yloc[YiU,1] ~ Yloc[YiU,2], col="green" )  # with covars and no other data issues
      points( Sloc[Si,1] ~ Sloc[Si,2], col="blue" ) # statistical locations
      points( p$plons[rcS[Si,1]] ~ p$plats[rcS[Si,2]] , col="purple", pch=25, cex=2 ) # check on rcS indexing
      points( p$plons[pa$Prow] ~ p$plats[ pa$Pcol] , col="cyan", pch=".", cex=0.01 ) # check on Proc Pcol indexing
      points( Ploc[pa$i,1] ~ Ploc[ pa$i, 2] , col="black", pch=20, cex=0.7 ) # check on pa$i indexing
    }

  # prediction stack:: check for covariates
    Pi = 1:nrow( Ploc ) 
    pa$i = Pi[ pa$i ]  ## flag to ignore
    kP = na.omit(pa$i)
    if ( length( kP) < p$n.min ) next()
    pa = pa[ which(is.finite( pa$i)),] 
    pa_n = nrow(pa)

    close(Ploc)

    # default output grid .. faaster form of expand.grid
    newdata = cbind( pa[ rep.int(1:pa_n, length(Ptime)), c("plon", "plat", "i")], rep.int(Ptime[], pa_n ))
    names(newdata) = c("plon", "plat", "i", "tiyr" )

    res = spacetime.timeseries ( p, YiU, Si, newdata )

    if ( is.null(res)) next()

    res = merge( res, pa, by="i", all.x=TRUE, all.y=FALSE )

    if (debugrun) {
      x11(); levelplot( mean ~ plon+plat, pa, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
      x11(); levelplot( sd   ~ plon+plat, pa, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }

    good = which( is.finite( rowSums(pa) ) )
    if (length(good) < 1) next()
    pa = pa[good,]

    # update P (predictions)
    ii = pa$i

    P = p$ff$P  # predictions
    test = rowSums( P[ii,] )
    u = which( is.finite( test ) )  # these have data already .. update
    if ( length( u ) > 0 ) {
      ui = ii[u]  # locations of P to modify
      Pn[ ui,,,] = Pn[ ui,,,] + 1 # update counts
      # update SD estimates of predictions with those from other locations via the
      # incremental  method ("online algorithm") of mean estimation after Knuth ;
      # see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
      stdev_update =  Psd[ ui,,, ] + ( pa$sd[u] -  Psd[ ui,,, ] ) / Pn[ui,,,]
      # update means: inverse-variance weighting   https://en.wikipedia.org/wiki/Inverse-variance_weighting
      means_update = ( P[ ui,,, ] / Psd[ ui,,, ]^2 + pa$mean[u] / pa$sd[u]^2 ) / ( Psd[ ui,,,]^(-2) + pa$sd[u]^(-2) )
      mm = which(is.finite( means_update + stdev_update ))
      if( length(mm)> 0) {
        # actual updates occur after everything has been computed first
        iumm = ui[mm]
        Psd[ iumm,,, ] = stdev_update[mm]
        P  [ iumm,,, ] = means_update[mm]
      }
    }

    # do this as a second pass in case NA's were introduced by the update .. unlikely , but just in case
    test = rowSums( P[ii,] )
    f = which( !is.finite( test ) ) # first time
    if ( length( f ) > 0 ) {
      fi = ii[f]
      Pn[ fi,,, ] = 1
      P[ fi,,, ] = pa$mean[f]
      Psd[ fi,,, ] = pa$sd[f]
    }

    #########

    if ( any( grepl ("statistics", p$spacetime.outputs))) {
      print( "Saving summary statisitics" )
      # save statistics last as this is an indicator of completion of all tasks .. restarts would be broken otherwise
      spacetime_stats = summarize model params
      S = p$ff$S  # statistical outputs inside loop to safely save data and pass onto other processes
      S[Si,1] = spacetime_stats["spatial error", "mode"]
      S[Si,2] = spacetime_stats["observation error", "mode"]
      S[Si,3] = spacetime_stats["range", "mode"]
      S[Si,4] = spacetime_stats["range", "sd"]
      if ( debugrun)  {
        print( spacetime_stats )
        cat( paste( Sys.time(), deid, "statistics saved  \n" ), file=p$debug.file, append=TRUE )
      }
      close(S)
    }

    if(0) {
      pps = expand.grid( plon=p$plons, plat=p$plats)
      # res = which(pps$plon > -50 & pps$plon < 50 & pps$plats < 50 & pps$plats > -50 ) # & P[,2] > 0   )
      res = which(pps$plon > min(pa$plon) & pps$plon < max(pa$plon) & pps$plat < max(pa$plat) & pps$plat > min(pa$plat) )
      x11(); levelplot( ( P[res,,,] ) ~ plon + plat, pps[res,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }

    rm( pa ) ; gc()
#    if ( debugrun) cat( paste( Sys.time(), deid, "end \n" ), file=p$debug.file, append=TRUE ) 
  }  # end for loop

  close(P)
  close(Y)

  return( "complete" )

}


