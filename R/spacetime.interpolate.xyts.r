

spacetime.interpolate.xyts = function( ip, p ) {

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  #---------------------
  # data for modelling
  # dependent vars # already link-transformed in spacetime.db("dependent")
  Y = p$ff$Y
  hasdata = 1:length(Y)
  bad = which( !is.finite( Y[]))
  if (length(bad)> 0 ) hasdata[bad] = NA

  # data locations
  Yloc = p$ff$Yloc
  bad = which( !is.finite( rowSums(Yloc[])))
  if (length(bad)> 0 ) hasdata[bad] = NA

  hasdata = na.omit(hasdata)

  #---------------------
  # prediction locations and covariates
  Ploc = p$ff$Ploc  # prediction locations
  phasdata = 1:nrow( Ploc ) # index of locs with no covariate data
  pbad = which( !is.finite( rowSums(Ploc[])))
  if (length(pbad)> 0 ) phasdata[ pbad ] = NA

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

  # main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    dd = p$runs[ iip, "jj" ]
    focal = t(Sloc[dd,])

    S = p$ff$S  # statistical outputs inside loop to safely save data and pass onto other processes
    if ( is.nan( S[dd,1] ) ) next()
    if ( !is.na( S[dd,1] ) ) next()
    S[dd,1] = NaN   # this is a flag such that if a run fails (e.g. in mesh generation), it does not get revisited
    # .. it gets over-written below if successful
    # choose a distance <= p$dist.max where n is within range of reasonable limits to permit a numerical solution
    # slow ... need to find a faster solution
    close(S)

    ppp = NULL
    ppp = try( point.in.block( focal[1,c(1,2)], Yloc[hasdata,], dist.max=p$dist.max, n.min=p$n.min, n.max=p$n.max,
      upsampling=p$upsampling, downsampling=p$downsampling, resize=TRUE ) )
    if( is.null(ppp)) next()
    if (class( ppp ) %in% "try-error" ) next()
    dist.cur = ppp$dist.to.nmax

    j = hasdata[ppp$indices]
    ndata = length(j) # number of data locations
    if (ndata < p$n.min) next()

    obs_ydata = list()
    obs_ydata[[ p$variables$Y ]] =  Y[j] # already link-transformed in spacetime.db("dependent")

    #  NOTE:: by default, all areas chosen to predict within the window.. but if covariates are involved,
    #  this can be done only where covariates exists .. so next step is to determine this and predict
    #  over the correct area.
    #  Once all predictions are complete, simple (kernal-based?) interpolation
    #  for areas without covariates can be completed
      windowsize.half = floor(dist.cur/p$pres)# convert distance to discretized increments of row/col indices
      pa_offsets = -windowsize.half : windowsize.half
      pa = expand.grid( Prow = rcS[dd,1] + pa_offsets, Pcol = rcS[dd,2] + pa_offsets ) # row,col coords
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
        plot( Yloc[ppp$indices,1]~ Yloc[ppp$indices,2], col="red", pch=".")
        points( Yloc[j,1] ~ Yloc[j,2], col="green" )
        points( focal[1,1] ~ focal[1,2], col="blue" )
        points( p$plons[rcS[dd,1]] ~ p$plats[rcS[dd,2]] , col="purple", pch=25, cex=2 )
        points( p$plons[pa$Prow] ~ p$plats[ pa$Pcol] , col="cyan", pch=".", cex=0.01 )
        points( Ploc[pa$i,1] ~ Ploc[ pa$i, 2] , col="yellow", pch=".", cex=0.7 )
      }

      close(Ploc)

      # prediction stack:: check for covariates
      if ( any( grepl ("predictions.direct", p$spacetime.outputs))) {
        pa$i = phasdata[ pa$i ]  ## flag to ignore
        kP = na.omit(pa$i)
        if ( length( kP) < p$n.min ) next()
        pa = pa[ which(is.finite( pa$i)),] # reduce to data locations as stack_index points only to locs with covariates
        preds_locs = as.matrix( pa[, c("plon", "plat")] )
        preds_eff = list()
        if ( exists( "COV", p$variables) ) {
          if ( length(p$variables$COV) == 1 ) {
            pcovars = as.data.frame(Pcov[ kP ])
          } else {
            pcovars = as.data.frame(Pcov[ kP,])
          }
          colnames( pcovars ) = p$variables$COV
          preds_eff[["covar"]] = as.list( pcovars )
        }
        preds_ydata = list()
        preds_ydata[[ p$variables$Y ]] = NA ## ie. to predict
      }


      # -----------
      # predictions
      if ( any( grepl ("predictions", p$spacetime.outputs))) {
        # do not use ifelse below ... it alters data structure
        print( "Prediction step ..")

        task = p$spacetime.outputs[grep("predictions", p$spacetime.outputs) ][1]  ## only take the first one in case there are many
        preds = NULL

      #  preds = .... model it here

        if (exists("spacetime.invlink", p)) {
          preds$xmean =  p$spacetime.invlink( preds$xmean )
          preds$xsd =  p$spacetime.invlink( preds$xsd )
        }

        if (is.null(preds)) {
          if ( debugrun) cat( paste( Sys.time(), deid, "prediction failure \n" ), file=p$debug.file, append=TRUE )
          next()
        }
        pa = cbind( pa, preds)
        if (debugrun) cat( paste(  Sys.time(), deid, "predictions completed \n" ), file=p$debug.file, append=TRUE )
        if (debugrun) {
          x11(); levelplot( xmean ~ plon+plat, pa, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
          x11(); levelplot( xsd   ~ plon+plat, pa, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        }
        good = which( is.finite( rowSums(pa) ) )
        if (length(good) < 1) next()
        pa = pa[good,]

        # update P (predictions)
        counts = 1 # column indices
        means = 2
        stdevs = 3
        ii = pa$i

        P = p$ff$P  # predictions
        test = rowSums( P[ii,] )
        u = which( is.finite( test ) )  # these have data already .. update
        if ( length( u ) > 0 ) {
          ui = ii[u]  # locations of P to modify
          P[ ui, counts ] = P[ ui, counts ] + 1 # update counts
          # update SD estimates of predictions with those from other locations via the
          # incremental  method ("online algorithm") of mean estimation after Knuth ;
          # see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
          stdev_update =  P[ ui, stdevs ] + ( pa$xsd[u] -  P[ ui, stdevs ] ) / P[ui, counts]
          # update means: inverse-variance weighting   https://en.wikipedia.org/wiki/Inverse-variance_weighting
          means_update = ( P[ ui, means ] / P[ ui, stdevs ]^2 + pa$xmean[u] / pa$xsd[u]^2 ) / ( P[ ui, stdevs]^(-2) + pa$xsd[u]^(-2) )
          mm = which(is.finite( means_update + stdev_update ))
          if( length(mm)> 0) {
            # actual updates occur after everything has been computed first
            P[ ui[mm], stdevs ] = stdev_update[mm]
            P[ ui[mm], means ]  = means_update[mm]
          }
        }

        # do this as a second pass in case NA's were introduced by the update .. unlikely , but just in case
        test = rowSums( P[ii,] )
        f = which( !is.finite( test ) ) # first time
        if ( length( f ) > 0 ) {
          fi = ii[f]
          P[ fi, counts ] = 1
          P[ fi, means ] = pa$xmean[f]
          P[ fi, stdevs ] = pa$xsd[f]
        }
      
      rm(MESH); gc()
    }

    #########

    if ( any( grepl ("statistics", p$spacetime.outputs))) {
      print( "Saving summary statisitics" )
      # save statistics last as this is an indicator of completion of all tasks .. restarts would be broken otherwise
      S = p$ff$S  # statistical outputs inside loop to safely save data and pass onto other processes
        S[dd,1] = inla.summary["spatial error", "mode"]
        S[dd,2] = inla.summary["observation error", "mode"]
        S[dd,3] = inla.summary["range", "mode"]
        S[dd,4] = inla.summary["range", "sd"]
        if ( debugrun)  {
          print( inla.summary )
          cat( paste( Sys.time(), deid, "statistics saved  \n" ), file=p$debug.file, append=TRUE )
        }
      close(S)
    }

    if(0) {
      pps = expand.grid( plon=p$plons, plat=p$plats)
      # zz = which(pps$plon > -50 & pps$plon < 50 & pps$plats < 50 & pps$plats > -50 ) # & P[,2] > 0   )
      zz = which(pps$plon > min(pa$plon) & pps$plon < max(pa$plon) & pps$plat < max(pa$plat) & pps$plat > min(pa$plat) )
      x11(); levelplot( ( P[zz,means] ) ~ plon + plat, pps[zz,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }


    rm( ii, good, pa, xs, xm, mi, mf, si, sf ) ; gc()
    if ( debugrun) cat( paste( Sys.time(), deid, "end \n" ), file=p$debug.file, append=TRUE ) 
  }  # end for loop

  close(P)
  close(Y)

  return( "complete" )

}


