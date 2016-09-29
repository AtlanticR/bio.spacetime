
  spacetime.interpolate.inla.local = function( ip=NULL, p, debugrun=FALSE ) {
    #\\ generic spatial and space-time interpolator using inla
    #\\ parameter and data requirements can be seen in bathymetry\src\bathymetry.r
    #\\ note this can run in parallel and serial mode

    # ip is the first parameter passed in the parallel mode
    if (exists( "libs", p)) RLibrary( p$libs )
    if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns
 
    # the following parameters are for inside and outside ... do not make them exact multiples as this seems to make things hang ..
    if ( !exists("inla.mesh.max.edge", stp))  stp$inla.mesh.max.edge = c(  0.025,   0.04 )    # proportion of 2*stp$dist.max or equivalent: c(inside,outside) -- must be positive valued
    if ( !exists("inla.mesh.offset", stp))  stp$inla.mesh.offset   = c( - 0.025,  - 0.05 )   # how much to extend inside and outside of boundary: proportion of dist.max .. neg val = proportion
    if ( !exists("inla.mesh.cutoff", stp)) stp$inla.mesh.cutoff   = c( - 0.05,   - 0.5 )    ## min distance allowed between points: proportion of dist.max ; neg val = proportion

    if ( !exists("inla.mesh.hull.radius", stp)) stp$inla.mesh.hull.radius = c( -0.04, - 0.08 ) ## radius of boundary finding algorythm ; neg val = proportion

    if ( !exists("inla.mesh.hull.resolution", stp)) stp$inla.mesh.hull.resolution = 125  ## resolution for discretization to find boundary

    if ( !exists("spacetime.noise", stp)) stp$spacetime.noise = 0.001  # add a little noise to coordinates to prevent a race condition

    if ( !exists("inla.alpha", stp)) stp$inla.alpha = 2 # bessel function curviness
    if ( !exists("inla.nsamples", stp)) stp$inla.nsamples = 5000 # posterior similations 
    if ( !exists("predict.in.one.go", stp)) stp$predict.in.one.go = FALSE # use false, one go is very very slow and a resource expensive method
    if ( !exists("predict.quantiles", stp)) stp$predict.quantiles = c(0.025, 0.975 )  # posterior predictions robustified by trimming extreme values 
    
    if ( !exists("debug.file", stp)) stp$debug.file = file.path( bio.workdirectory, "inla.debug.out" )

    # priors
    # kappa0 = sqrt(8)/p$expected.range
    # tau0 = 1/(sqrt(4*pi)*kappa0* p$expected.sigma)

    #---------------------
    # data for modelling
    # dependent vars # already link-transformed in spacetime.db("dependent")
    Y =  p$ff$Y  # readonly
    Yi = 1:length(Y)
    bad = which( !is.finite( Y[]))
    if (length(bad)> 0 ) Yi[bad] = NA

   # data locations
    Yloc =  p$ff$Yloc
    bad = which( !is.finite( rowSums(Yloc[])))
    if (length(bad)> 0 ) Yi[bad] = NA

    # covariates (independent vars)
    if ( exists( "COV", p$variables) ) {
      Ycov =  p$ff$Ycov
      if ( length( p$variables$COV ) == 1 ) {
        bad = which( !is.finite( Ycov[]) )
      } else {
        bad = which( !is.finite( rowSums(Ycov[])) )
      }
      if (length(bad)> 0 ) Yi[bad] = NA
    }

    Yi = na.omit( Yi )

    #---------------------
    # prediction locations and covariates
    Ploc =  p$ff$Ploc # read only
    Pi = 1:nrow( Ploc ) # index of locs with no covariate data
    pbad = which( !is.finite( rowSums(Ploc[])))
    if (length(pbad)> 0 ) Pi[ pbad ] = NA
    if ( exists( "COV", p$variables) ) {
      Pcov =  p$ff$Pcov  # covariates at prediction locations; read only
      if ( length( p$variables$COV ) == 1 ) {
        pbad = which( !is.finite( Pcov[]) )
      } else {
        pbad = which( !is.finite( rowSums(Pcov[])) )
      }
      if (length(pbad)> 0 ) Pi[pbad] = NA
      close(Pcov)
    }
    rcP = data.frame( cbind( Prow = (Ploc[,1]-p$plons[1])/p$pres + 1,  Pcol = (Ploc[,2]-p$plats[1])/p$pres + 1))
    # rcP$i =1:nrow(rcP) # row index
    rcP$rc = paste( rcP$Prow, rcP$Pcol, sep="~")
    rcP$Prow = rcP$Pcol = NULL
    gc()

    #-----------------
    # row, col indices
    Sloc =  p$ff$Sloc  # statistical output locations
    rcS = data.frame( cbind( Srow = (Sloc[,1]-p$plons[1])/p$pres + 1,  Scol = (Sloc[,2]-p$plats[1])/p$pres + 1))

    # main loop over each output location in S (stats output locations)
    for ( iip in ip ) {
      Si = p$runs[ iip, "locs" ]
      # Si=2000
      if (debugrun) deid = paste( Sys.info()["nodename"], "index=", Si )
      if (debugrun) cat( paste( Sys.time(), deid, "start \n" ), file=p$debug.file, append=TRUE )
      focal = t(Sloc[Si,])

      S = p$ff$S  # statistical outputs
      if ( is.nan( S[Si,1] ) ) next()
      if ( !is.na( S[Si,1] ) ) next()
      S[Si,1] = NaN   # this is a flag such that if a run fails (e.g. in mesh generation), it does not get revisited
      # .. it gets over-written below if successful
      # choose a distance <= p$dist.max where n is within range of reasonable limits to permit a numerical solution
      close(S)

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

      if (debugrun) cat( paste( Sys.time(), deid, "n=", ndata, "dist=", dist.cur, "\n" ), file=p$debug.file, append=TRUE )
      locs_noise = runif( ndata*2, min=-p$pres*p$spacetime.noise, max=p$pres*p$spacetime.noise ) # add  noise  to prevent a race condition

      # also sending direct distances rather than proportion seems to cause issues..
      MESH = spacetime.mesh( locs=Yloc[YiU,]+locs_noise,
        lengthscale=dist.cur*2,
        max.edge=p$inla.mesh.max.edge * dist.cur*2,
        bnd.offset=p$inla.mesh.offset,
        cutoff=p$inla.mesh.cutoff,
        convex=p$inla.mesh.hull.radius,
        resolution=p$inla.mesh.hull.resolution )

      rm(locs_noise)

      if ( debugrun)  {
        if ( is.null(MESH)) {
          cat( paste( Sys.time(), deid,  "mesh try error \n" ), file=p$debug.file, append=TRUE )
          next()
        } else {
          cat( paste( Sys.time(), deid,  "mesh finished \n" ), file=p$debug.file, append=TRUE )
          plot( MESH )  # or ... spacetime.plot( p=p, "mesh", MESH=MESH )
        }
      }


      SPDE = inla.spde2.matern( MESH,
        alpha=p$inla.alpha #, # alpha is the Bessel smoothness factor .. 1(?) gives exponential correlation function
        # B.tau=matrix(c(log(tau0),-1,+1),nrow=1,ncol=3),
        # B.kappa=matrix(c(log(kappa0),0,-1),nrow=1,ncol=3),
        # theta.prior.mean=c(0,0), # thetas are the internal representations of prior offsets for tau and kappa (i.e.,range and sigma)
        # theta.prior.prec=c(0.1, 0.1) # precision of thetas
      )

      # effects .. a list of two elements: first is for SPDE and second is for covariates
      obs_index = inla.spde.make.index(name=p$spatial.field.name, SPDE$n.spde)
      obs_eff = list()
      obs_eff[["spde"]] = c( obs_index, list(intercept=1) )
      if ( exists( "COV", p$variables) ) {
        if ( length(p$variables$COV) == 1 ) {
          covar = as.data.frame( Ycov[YiU])
        } else {
          covar= as.data.frame( Ycov[YiU,])
        }
        colnames( covar ) = p$variables$COV
        obs_eff[["covar"]] = as.list(covar)
        obs_A = list( inla.spde.make.A( mesh=MESH, loc=Yloc[YiU,] ), 1 )
      } else {
        obs_A = list( inla.spde.make.A( mesh=MESH, loc=Yloc[YiU,] ) ) # no effects
      }

      obs_ydata = list()
      obs_ydata[[ p$variables$Y ]] = Y[YiU]
      if ( exists("spacetime.link", p) ) obs_ydata[[ p$variables$Y ]] = p$spacetime.link ( Y[YiU] ) 
      
      DATA = inla.stack( tag="obs", data=obs_ydata, A=obs_A, effects=obs_eff, remove.unused=FALSE )
      rm ( obs_index, obs_eff, obs_ydata, obs_A )
      # remove.unused=FALSE ensures that we can look at the estimated field effect without
      #   having to do expensive separate predictions.
      # DATA$A is projection matrix to translate from mesh nodes to data nodes

      # -----------
      # PREDICTIONS
      #      NOTE:: by default, all areas chosen to predict within the window.. but if covariates are involved,
      #        this can be done only where covariates exists .. so next step is to determine this and predict
      #        over the correct area.
      #        Once all predictions are complete, simple (kernal-based?) interpolation
      #        for areas without covariates can be completed
      windowsize.half = floor(dist.cur/p$pres)# covert distance to discretized increments of row/col indices
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
        plot( Yloc[U,1]~ Yloc[U,2], col="red", pch=".")
        points( Yloc[YiU,1] ~ Yloc[YiU,2], col="green" )
        points( focal[1,1] ~ focal[1,2], col="blue" )
        points( p$plons[rcS[Si,1]] ~ p$plats[rcS[Si,2]] , col="purple", pch=25, cex=2 )
        points( p$plons[pa$Prow] ~ p$plats[ pa$Pcol] , col="cyan", pch=".", cex=0.01 )
        points( Ploc[pa$i,1] ~ Ploc[ pa$i, 2] , col="yellow", pch=".", cex=0.7 )
      }

      # prediction stack:: check for covariates
      if ( any( grepl ("predictions.direct", p$spacetime.outputs))) {
        pa$i = Pi[ pa$i ]  ## flag to ignore
        kP = na.omit(pa$i)
        if ( length( kP) < p$n.min ) next()
        pa = pa[ which(is.finite( pa$i)),] # reduce to data locations as stack_index points only to locs with covariates
        preds_locs = as.matrix( pa[, c("plon", "plat")] )
        preds_index = inla.spde.make.index( name=p$spatial.field.name, SPDE$n.spde)
        preds_eff = list()
        preds_eff[["spde"]] = c( preds_index, list(intercept=1) )
        if ( exists( "COV", p$variables) ) {
          Pcov =  p$ff$Pcov  # covariates at prediction locations; read only
          if ( length(p$variables$COV) == 1 ) {
            pcovars = as.data.frame(Pcov[ kP ])
          } else {
            pcovars = as.data.frame(Pcov[ kP,])
          }
          colnames( pcovars ) = p$variables$COV
          preds_eff[["covar"]] = as.list( pcovars )
          preds_A = list( inla.spde.make.A(MESH, loc=preds_locs ), 1)
          close(Pcov)
        } else {
          preds_A = list( inla.spde.make.A(MESH, loc=preds_locs ) )
        }
        preds_ydata = list()
        preds_ydata[[ p$variables$Y ]] = NA ## ie. to predict
        PREDS = inla.stack( tag="preds", data=preds_ydata, A=preds_A, effects=preds_eff, remove.unused=FALSE )
        DATA = inla.stack(DATA, PREDS )
        preds_stack_index = inla.stack.index( DATA, "preds")$data  # indices of predictions in stacked data
        rm ( preds_eff, preds_ydata, preds_A, PREDS, preds_index, preds_locs, pcovars, kP ); gc()
      }

      close(Yloc)
      close(Ploc)

      RES = NULL
      RES = spacetime.inla.call( FM=p$modelformula, DATA=DATA, SPDE=SPDE, FAMILY=p$spacetime.family )
      if ( debugrun && !is.null( RES) ) {
        cat( paste(  Sys.time(), deid, "computations finished \n" ), file=p$debug.file, append=TRUE )
        print(RES)
        print( summary(RES))
        # low level debugging .. and looking at posterior marginals
        idat =  inla.stack.index( DATA, 'data')$data # indices of data locations
        spacetime.plot( p=p, "range", RES=RES, MESH=MESH, SPDE=SPDE, vname=p$spatial.field.name, idat=idat )
        spacetime.plot( p=p, "nugget", RES=RES, MESH=MESH, SPDE=SPDE, vname=p$spatial.field.name, idat=idat )
        spacetime.plot( p=p, "partial.sill", RES=RES, MESH=MESH, SPDE=SPDE, vname=p$spatial.field.name, idat=idat )
        spacetime.plot( p=p, "fixed.intercept", RES=RES, MESH=MESH, SPDE=SPDE, vname=p$spatial.field.name, idat=idat )
        rm (idat); gc()
      }

      if (is.null(RES)) {
        cat( paste(  Sys.time(), Sys.info()["nodename"], "index=", Si,  "inla call error \n" ), file=p$debug.file, append=TRUE )
        next()
      }

      # inla.spde2.matern creates files to disk that are not cleaned up:
      spdetmpfn = SPDE$f$spde2.prefix
      fns = list.files( dirname( spdetmpfn ), all.files=TRUE, full.names=TRUE, recursive=TRUE, include.dirs=TRUE )
      oo = grep( basename(spdetmpfn), fns )
      if(length(oo)>0) file.remove( sort(fns[oo], decreasing=TRUE) )
      rm ( DATA, oo, fns, spdetmpfn); gc()

      # -----------
      # predictions
      if ( any( grepl ("predictions", p$spacetime.outputs))) {
        # do not use ifelse below ... it alters data structure
        print( "Prediction step ..")

        task = p$spacetime.outputs[grep("predictions", p$spacetime.outputs) ][1]  ## only take the first one in case there are many
        preds = NULL

        if ( task=="predictions.direct" ) {
          # precomputed ... slow and expensive in RAM/CPU, just extract from tag indices
          xmean = RES$summary.fitted.values[ preds_stack_index, "mean"]
          xsd = RES$summary.fitted.values[ preds_stack_index, "sd"]
          if (exists("spacetime.invlink", p)) {
            xmean =  p$spacetime.invlink( xmean )
            xsd =  p$spacetime.invlink( xsd )
          }
          preds = data.frame(xmean=xmean, xsd=xsd) 
        }

        if ( task=="predictions.projected") {
          #\\ note this method only works with simple aSiitive models 
          #\\ when smoothes are involved, it becomes very complicated and direct estimation is probably faster/easier
          locs_new=pa[,c("plon", "plat")]
          pG = inla.mesh.projector( MESH, loc=as.matrix( locs_new) )
          posterior.samples = inla.posterior.sample(n=p$inla.nsamples, RES)
          rnm = rownames(posterior.samples[[1]]$latent )  
          posterior = sapply( posterior.samples, p$spacetime.posterior.extract, rnm=rnm )
          if (exists("spacetime.invlink", p)) posterior = p$spacetime.invlink( posterior )   # return to user scale
          rm(posterior.samples); gc()
          # robustify the predictions by trimming extreme values .. will have minimal effect upon mean
          # but variance estimates should be useful/more stable as the tails are sometimes quite long 
          for (ii in 1:nrow(posterior )) {
            qnt = quantile( posterior[ii,], probs=p$predict.quantiles, na.rm=TRUE ) 
            toolow = which( posterior[ii,] < qnt[1] )
            toohigh = which (posterior[ii,] > qnt[2] )
            if (length( toolow) > 0 ) posterior[ii,toolow] = qnt[1]
            if (length( toohigh) > 0 ) posterior[ii,toohigh] = qnt[2]
          }
          xmean = c( inla.mesh.project( pG, field=apply( posterior, 1, mean, na.rm=TRUE )  ))
          xsd   = c( inla.mesh.project( pG, field=apply( posterior, 1, sd, na.rm=TRUE )  ))
          preds = data.frame(xmean=xmean, xsd=xsd) 
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

        P =  p$ff$P # predictions
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
        close(P)

      }
      rm(MESH); gc()

      # bathymetry.figures( DS="predictions", p=p ) # to debug

      if ( any( grepl ("statistics", p$spacetime.outputs))) {
        print( "Saving summary statisitics" )
        # extract summary statistics from a spatial (SPDE) analysis and update the output file
        inla.summary = spacetime.summary.inla.spde2 ( RES, SPDE )
        # save statistics last as this is an indicator of completion of all tasks .. restarts would be broken otherwise

        S = p$ff$S  # statistical outputs
        S[Si,1] = inla.summary["spatial error", "mode"]
        S[Si,2] = inla.summary["observation error", "mode"]
        S[Si,3] = inla.summary["range", "mode"]
        S[Si,4] = inla.summary["range", "sd"]
        if ( debugrun)  {
          print( inla.summary )
          cat( paste( Sys.time(), deid, "statistics saved  \n" ), file=p$debug.file, append=TRUE )
        }
        close(S)
      }

      if(debugrun) {
        P = p$ff$P
        pps = expand.grid( plon=p$plons, plat=p$plats)
        # zz = which(pps$plon > -50 & pps$plon < 50 & pps$plats < 50 & pps$plats > -50 ) # & P[,2] > 0   )
        zz = which(pps$plon > min(pa$plon) & pps$plon < max(pa$plon) & pps$plat < max(pa$plat) & pps$plat > min(pa$plat) )
        x11(); levelplot( ( P[zz,means] ) ~ plon + plat, pps[zz,], aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        close(P)
      }

      close(Y)

      rm( SPDE, RES) ; gc()

      rm( ii, good, pa, xs, xm, mi, mf, si, sf ); gc()
      if ( debugrun) cat( paste( Sys.time(), deid, "end \n" ), file=p$debug.file, append=TRUE )
    
    }  # end for loop
    return( "complete" )
  }


