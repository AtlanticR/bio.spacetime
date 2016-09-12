

spacetime.interpolate.gam.harmonic = function( ip, p ) {
  #\\ harmonic in time method

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns
   # load hdf5 data objects pointers
  p = spacetime.db( p=p, DS="filenames" )



if(0) {

  mf = switch( p$tsmethod ,
    annual = ' t ~ s(yr) ',
    seasonal.basic = ' t ~ s(yr) + s(dyear, bs="cc") ',
    seasonal.smoothed = ' t ~ s(yr, dyear) + s(yr) + s(dyear, bs="cc")  ',
    harmonics.1 = ' t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w)  ',
    harmonics.2 = ' t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) ' ,
    harmonics.3 = ' t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w) + s(yr, cos.w2) + s(yr, sin.w2) + s(cos.w2) + s( sin.w2 ) + s(yr, cos.w3) + s(yr, sin.w3)  + s(cos.w3) + s( sin.w3 ) '
  )

  mf = formula(mf)

  for ( dm in p$dist.km ) {
      drange = c(-1,1) * dm
      plon0 = pp$plon + drange
      plat0 = pp$plat + drange
      i = which( bb$plon > plon0[1] & bb$plon < plon0[2] & bb$plat > plat0[1] & bb$plat < plat0[2] )
      if (length(i) > p$nMin.tbot ) {

      #  browser()

        # only attempt interpolation if we have enough data (nMin.tbot)
        x = bb[i,] # faster to reduce the size of bb here
        # remove potentially noisy/erroneous data --- they are highly influential when there is little data
        # xt = quantile( x$t, probs=c(0.005, 0.995) )
        # xi = which( x$t >= xt[1] & x$t <= xt[2] )
        # if (length(xi) < p$nMin.tbot ) next()
        # x = x[xi, ]

        x$w = 1 / (( pp$plon - x$plon)**2 + (pp$plat - x$plat)**2 )# weight data in space: inverse distance squared
        x$w[ which( x$w < 1e-3 ) ] = 1e-3
        x$w[ which( x$w > 1 ) ] = 1

        # data transformations and creation of new variables where required for raw data
        if ( p$tsmethod %in% c( "harmonics.1", "harmonics.2", "harmonics.3"  ) ) {
          x$cos.w  = cos( 2*pi*x$tiyr )
          x$sin.w  = sin( 2*pi*x$tiyr )

          years.with.data = unique( x$yr)
          no.years = which( !( zz$yr %in% years.with.data) )
          zz$yr[ no.years ] = median( years.with.data)
          zz$cos.w  = cos( zz$tiyr )
          zz$sin.w  = sin( zz$tiyr )

          # compute additional harmonics only if required (to try to speed things up a bit)
          if ( p$tsmethod %in% c( "harmonics.2", "harmonics.3"  ) ) {
            x$cos.w2 = cos( 2*x$tiyr )
            x$sin.w2 = sin( 2*x$tiyr )
            zz$cos.w2 = cos( 2*zz$tiyr )
            zz$sin.w2 = sin( 2*zz$tiyr )
          }
          if ( p$tsmethod %in% c( "harmonics.3"  ) ) {
            x$cos.w3 = cos( 3*x$tiyr )
            x$sin.w3 = sin( 3*x$tiyr )
            zz$cos.w3 = cos( 3*zz$tiyr )
            zz$sin.w3 = sin( 3*zz$tiyr )
          }
        }

        tsmodel = NULL
        tsmodel = try( gam( mf, data=x, weights=w, optimizer=c("outer","bfgs")  ) )
        if ( ! "try-error" %in% class(tsmodel) ) {
          out = try( predict( tsmodel, newdata=zz, type="response", se.fit=T ) )
          if ( ! "try-error" %in% class( out ) ) {
            zz$fit = out$fit
            zz$se = out$se.fit
            break()  # candidate predictions found exit inner loop (dm)
          }
        }
      } 
  }
}


  #---------------------
  # data for modelling
  # dependent vars # already link-transformed in spacetime.db("dependent")
  Y = h5file(p$ptr$Y)["Y"]
  hasdata = 1:length(Y)
  bad = which( !is.finite( Y[]))
  if (length(bad)> 0 ) hasdata[bad] = NA

  # data locations
  Yloc = h5file(p$ptr$Yloc)["Yloc"]
  bad = which( !is.finite( rowSums(Yloc[])))
  if (length(bad)> 0 ) hasdata[bad] = NA

  # covariates (independent vars)
  if ( exists( "COV", p$variables) ) {
    Ycov = h5file(p$ptr$Ycov)["Ycov"]
    if ( length( p$variables$COV ) == 1 ) {
      bad = which( !is.finite( Ycov[]) )
    } else {
      bad = which( !is.finite( rowSums(Ycov[])) )
    }
    if (length(bad)> 0 ) hasdata[bad] = NA
  }

  hasdata = na.omit(hasdata)

  #---------------------
  # prediction locations and covariates
  Ploc = h5file(p$ptr$Ploc)["Ploc"]  # prediction locations
  phasdata = 1:nrow( Ploc ) # index of locs with no covariate data
  pbad = which( !is.finite( rowSums(Ploc[])))
  if (length(pbad)> 0 ) phasdata[ pbad ] = NA
  if ( exists( "COV", p$variables) ) {
    Pcov = h5file(p$ptr$Pcov)["Pcov"]  # covariates at prediction locations
    if ( length( p$variables$COV ) == 1 ) {
      pbad = which( !is.finite( Pcov[]) )
    } else {
      pbad = which( !is.finite( rowSums(Pcov[])) )
    }
    if (length(pbad)> 0 ) phasdata[pbad] = NA
  }
  rcP = data.frame( cbind( Prow = (Ploc[,1]-p$plons[1])/p$pres + 1,  Pcol = (Ploc[,2]-p$plats[1])/p$pres + 1))
  # rcP$i =1:nrow(rcP) # row index
  rcP$rc = paste( rcP$Prow, rcP$Pcol, sep="~")
  rcP$Prow = rcP$Pcol = NULL
  gc()

  #-----------------
  # row, col indices
  Sloc = h5file(p$ptr$Sloc)["Sloc"]  # statistical output locations
  rcS = data.frame( cbind( Srow = (Sloc[,1]-p$plons[1])/p$pres + 1,  Scol = (Sloc[,2]-p$plats[1])/p$pres + 1))

  # main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    dd = p$runs[ iip, "jj" ]
    focal = t(Sloc[dd,])

    S = h5file(p$ptr$S)["S"]  # statistical outputs

    if ( is.nan( S[dd,1] ) ) next()
    if ( !is.na( S[dd,1] ) ) next()

    S[dd,1] = NaN   # this is a flag such that if a run fails (e.g. in mesh generation), it does not get revisited
    # .. it gets over-written below if successful
    # choose a distance <= p$dist.max where n is within range of reasonable limits to permit a numerical solution
    # slow ... need to find a faster solution
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

        P = h5file(p$ptr$P)["P"]  # predictions
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
      S[dd,1] = inla.summary["range", "mode"]
      S[dd,2] = inla.summary["range", "sd"]
      S[dd,3] = inla.summary["spatial error", "mode"]
      S[dd,4] = inla.summary["observation error", "mode"]
      if ( debugrun)  {
        print( inla.summary )
        cat( paste( Sys.time(), deid, "statistics saved  \n" ), file=p$debug.file, append=TRUE )
      }
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

  return( "complete" )

}


