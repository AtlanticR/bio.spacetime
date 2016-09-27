

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

  # data locations
  Ytime = p$ff$Ytime
  bad = which( !is.finite( rowSums(Ytime[])))
  if (length(bad)> 0 ) hasdata[bad] = NA

  hasdata = na.omit(hasdata)

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
    dd = p$runs[ iip, "locs" ]
    focal = Sloc[dd,]

    S = p$ff$S  # statistical outputs inside loop to safely save data and pass onto other processes
    if ( is.nan( S[dd,1] ) ) next()
    if ( !is.na( S[dd,1] ) ) next()
    S[dd,1] = NaN   # this is a flag such that if a run fails (e.g. in mesh generation), it does not get revisited
    # .. it gets over-written below if successful
    # choose a distance <= p$p$dist.max where n is within range of reasonable limits to permit a numerical solution
    # slow ... need to find a faster solution
    close(S)

    dlon = abs( Sloc[dd,1] - Yloc[hasdata,1] ) 
    dlat = abs( Sloc[dd,2] - Yloc[hasdata,2] ) 

    U = NULL  # data within a specified range/distance/number of points
    U = which( dlon  <= p$dist.max  & dlat <= p$dist.max ) # faster to take a block 
    ndat = length(U)
    dist.cur = p$dist.max

    # moved point in polygon here as operating on ff means simpler usage of pointers
    # rather than working with a copy of the data    
    if ( ndat < p$n.min )  {
      for ( usamp in p$upsampling )  {
        dist.cur = p$dist.max * usamp
        U = which( dlon < dist.cur & dlat < dist.cur ) # faster to take a block 
        ndat = length(U)
        if ( ndat >= p$n.min ) {
          if (ndat >= p$n.max) {
            U = U[ .Internal( sample(  length(U), p$n.max, replace=FALSE, prob=NULL)) ] 
          } 
          break()
        }
      }
    } else if ( ndat >= p$n.min ) {
      if ( ndat <= p$n.max * 1.5 ) { # if close to p$n.max, subsample quickly and return
        if ( ndat > p$n.max) { 
          U = U[ .Internal( sample(  length(U), p$n.max, replace=FALSE, prob=NULL)) ] 
          break()
        }        
      } else {
        for ( dsamp in p$downsampling )  {
          dist.cur = p$dist.max * dsamp
          U = which( dlon < dist.cur & dlat < dist.cur )# faster to take a block 
          ndat = length(U)
          if ( ndat <= p$n.max ) break()
          if ( dist.cur <= dist.min ) {
            # reached lower limit in distance, taking a subsample instead
            U = which( dlon < dist.min & dlat < dist.min ) # faster to take a block 
            U = U[ .Internal( sample( length(U), p$n.max, replace=FALSE, prob=NULL)) ]
            break()
          }
        }
        # as a last resort, try sampling from the data as the distances are getting too small
        U = U[ .Internal( sample( length(U), p$n.max, replace=FALSE, prob=NULL)) ]
        if ( ndat <= p$n.max ) break()
      }
    }

    if ( is.null(U)) next()
    ndat = length(U)
    if ((ndat < p$n.min) | (ndat > p$n.max) ) next()

    j = hasdata[U]
    ndata = length(j) # number of data locations
    if (ndata < p$n.min) next()

# -->> p$nMin.tbot == p$n.min

    #  NOTE:: by default, all areas chosen to predict within the window.. but if covariates are involved,
    #  this can be done only where covariates exists .. so next step is to determine this and predict
    #  over the correct area.
    #  Once all predictions are complete, simple (kernal-based?) interpolation
    #  for areas without covariates can be completed
      windowsize.half = floor(dist.cur/p$pres) # convert distance to discretized increments of row/col indices
      pa_offsets = -windowsize.half : windowsize.half
      pa = expand.grid( Prow = rcS[dd,1] + pa_offsets, Pcol = rcS[dd,2] + pa_offsets ) # row,col coords
   
      pa$Prow = rcS[dd,1] + pa_offsets
      pa$Pcol = rcS[dd,2] + pa_offsets  # row,col coords
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
        plot( Yloc[U,1]~ Yloc[U,2], col="red", pch=".")
        points( Yloc[j,1] ~ Yloc[j,2], col="green" )
        points( Sloc[dd,1] ~ Sloc[dd,2], col="blue" )
        points( p$plons[rcS[dd,1]] ~ p$plats[rcS[dd,2]] , col="purple", pch=25, cex=2 )
        points( p$plons[pa$Prow] ~ p$plats[ pa$Pcol] , col="cyan", pch=".", cex=0.01 )
        points( Ploc[pa$i,1] ~ Ploc[ pa$i, 2] , col="yellow", pch=".", cex=0.7 )
      }

      close(Ploc)

    # prediction stack:: check for covariates
      phasdata = 1:nrow( Ploc )
      pa$i = phasdata[ pa$i ]  ## flag to ignore
      kP = na.omit(pa$i)
      if ( length( kP) < p$n.min ) next()
      pa = pa[ which(is.finite( pa$i)),] # reduce to data locations as stack_index points only to locs with covariates

              # default output grid
 
              res = expand.grid( Ptime , Ploc )

              # harmonic method explanation:
              # old method use 1 harmonic ... forcing sinusoid as a seasonal component
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
            

                x = Y[j,] # faster to reduce the size of Y here
                # ---> Y[j] # already link-transformed in spacetime.db("dependent")
                
                if(0) {
                  # is using inverse distance weighting
                  Y_wgt = 1 / (( Sloc[dd,1] - Yloc[j,1])**2 + (Sloc[dd,2] - Yloc[j,2])**2 )# weight data in space: inverse distance squared
                  Y_wgt[ which( Y_wgt < 1e-3 ) ] = 1e-3
                  Y_wgt[ which( Y_wgt > 1 ) ] = 1
                }
               
                # data transformations and creation of new variables where required for raw data 
                if ( p$tsmethod %in% c( "harmonics.1", "harmonics.2", "harmonics.3"  ) ) {
                  Y_cos_w  = cos( 2*pi*x$tiyr )
                  Y_sin_w  = sin( 2*pi*x$tiyr )
                 
                  years.with.data = unique( x$yr)
                  no.years = which( !( res$yr %in% years.with.data) )
                  res$yr[ no.years ] = median( years.with.data) 
                  res$cos.w  = cos( res$tiyr )
                  res$sin.w  = sin( res$tiyr )
                  
                  # compute additional harmonics only if required (to try to speed things up a bit)
                  if ( p$tsmethod %in% c( "harmonics.2", "harmonics.3"  ) ) {
                    Y_cos_w2 = cos( 2*x$tiyr )
                    Y_sin_w2 = sin( 2*x$tiyr )
                    res$cos.w2 = cos( 2*res$tiyr )
                    res$sin.w2 = sin( 2*res$tiyr )
                  }
                  if ( p$tsmethod %in% c( "harmonics.3"  ) ) {
                    Y_cos_w3 = cos( 3*x$tiyr )
                    Y_sin_w3 = sin( 3*x$tiyr )
                    res$cos.w3 = cos( 3*res$tiyr )
                    res$sin.w3 = sin( 3*res$tiyr )
                  }
                }

                # estimate model parameters
                tsmodel = try( 
                  gam( t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w), 
                    data=x, weights=Y_wgt, optimizer=c("outer","bfgs")  ) ) 
              
                if ( class(tsmodel) %in% "try-error" ) next()
               
                out = try( predict( tsmodel, newdata=res, type="response", se.fit=T ) ) 
                  if ( ! "try-error" %in% class( out ) ) {
                    res$fit = out$fit
                    res$se = out$se.fit 
                  }
                }
              }
               if ( any(is.finite(res$fit)) ) {
                print (iip)
                P[ iip,] = res$fit
                P.se[iip,] = res$se
              }
            } # end each point


      if (is.null(preds)) {
        if ( debugrun) cat( paste( Sys.time(), deid, "prediction failure \n" ), file=p$debug.file, append=TRUE )
        next()
      }
      if (exists("spacetime.invlink", p)) {
        preds$xmean =  p$spacetime.invlink( preds$xmean )
        preds$xsd =  p$spacetime.invlink( preds$xsd )
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
      # res = which(pps$plon > -50 & pps$plon < 50 & pps$plats < 50 & pps$plats > -50 ) # & P[,2] > 0   )
      res = which(pps$plon > min(pa$plon) & pps$plon < max(pa$plon) & pps$plat < max(pa$plat) & pps$plat > min(pa$plat) )
      x11(); levelplot( ( P[res,means] ) ~ plon + plat, pps[res,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }

    rm( ii, good, pa, xs, xm, mi, mf, si, sf ) ; gc()
    if ( debugrun) cat( paste( Sys.time(), deid, "end \n" ), file=p$debug.file, append=TRUE ) 
  }  # end for loop

  close(P)
  close(Y)

  return( "complete" )

}


