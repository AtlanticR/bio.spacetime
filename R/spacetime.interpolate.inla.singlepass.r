
  spacetime.interpolate.inla.singlepass = function( FM="Y ~ -1 + intercept + f( spatial.field, model=SPDE )",
    Y, locs, plocs, covar=NULL, pcovar=NULL, method="fast", link="identity", nsamples=5000 ) {
    #\\ low-level function -- single pass/fast spatial interpolator using inla
    #\\   but no bigmemory objects and no stats, etc. .. prediction only, serial mode only
    #\\   method="fast" using an indirect estimate based upon posterior projections of the input
    #\\   method="direct" uses direct (more accurate) estimation but way too slow and resource
    #\\     intensive to be useful for production runs .. used only to check the fast method
    #\\ default formula is a simple RF no covariates

    require(INLA)

    if (0) {
       # just for debugging / testing ... and example of access method:
      bioLibrary( "bio.utilities", "bio.spacetime")
        require(sp)
        data(meuse)
        x = meuse[, c("x")]
        y = meuse[, c("y")]
        z = meuse$zinc
        plocs = expand.grid(
          seq( min(x), max(x), length.out=100),
          seq( min(y), max(y), length.out=100), KEEP.OUT.ATTRS=FALSE )
        names(plocs) = colnames(locs)

        covar=NULL
        method="fast"
        link="log"
        link ="identity"
        nsamples=1000
        FM="Y ~ -1 + intercept + f( spatial.field, model=SPDE )"
        res = spacetime.interpolate.inla.singlepass (
          Y=z, locs=cbind(x,y), plocs=plocs, method=method, link=link )

        require(lattice)
        levelplot( xmean ~ x+y, data=res )

    }

    # identity links by default .. add more if needed here
    spacetime.link = function(X) {X}
    spacetime.invlink = function(X) {X}

    if (link=="log" ) {
      spacetime.link = function(X) {log(X)}
      spacetime.invlink = function(X) {exp(X)}
    }

    locs = as.matrix( locs)

    lengthscale = max( diff(range( locs[,1])), diff(range( locs[,2]) )) / 10  # in absence of range estimate take 1/10 of domain size

    ndata = length(Y)
    noise = lengthscale * 1e-9
    locs = locs + runif( ndata*2, min=-noise, max=noise ) # add  noise  to prevent a race condition .. inla does not like uniform grids

    MESH = spacetime.mesh( locs, lengthscale=lengthscale )
    if ( is.null( MESH) ) return( "Mesh Error" )

    SPDE = inla.spde2.matern( MESH,  alpha=2 ) # alpha is 2*nu (Bessel smoothness factor)

    # FM = formula( Y ~ -1 + intercept + f( spatial.field, model=SPDE ) )
    FM = formula( FM )
    varY = as.character( FM[2] )
    obs_index = inla.spde.make.index(name="spatial.field", SPDE$n.spde)
    obs_eff = list()
    obs_eff[["spde"]] = c( obs_index, list(intercept=1) )

    if ( !is.null(covar) ) {
        covar = as.data.frame( covar[,])
        obs_eff[["covar"]] = as.list(covar)
        obs_A = list( inla.spde.make.A( mesh=MESH, loc=locs[,] ), 1 )
    } else {
        obs_A = list( inla.spde.make.A( mesh=MESH, loc=locs[,] ) ) # no effects
    }

    obs_ydata = list()
    obs_ydata[[ varY ]] = spacetime.link ( Y )
    DATA = inla.stack( tag="obs", data=obs_ydata, A=obs_A, effects=obs_eff, remove.unused=FALSE )

    if ( method=="direct") {
      # direct method
      preds_index = inla.spde.make.index( name="spatial.field", SPDE$n.spde)
      preds_eff = list()
      preds_eff[["spde"]] = c( list( intercept=rep(1,MESH$n )),
           inla.spde.make.index(name="spatial.field", n.spde=SPDE$n.spde) )
      if ( !is.null( pcovar) ) preds_eff[["covar"]] = as.list( as.data.frame( pcovar ))
      ydata = list()
      ydata[[ varY ]] = NA
      Apreds = inla.spde.make.A(MESH, loc=as.matrix( plocs ) )
      PREDS = inla.stack( tag="preds",
        data=list( Y=NA),
        A=list(Apreds),
        effects=list(
          c( list(intercept=rep(1, MESH$n)),
             inla.spde.make.index( name="spatial.field", MESH$n)))
      )
      DATA = inla.stack(DATA, PREDS)
      i_data = inla.stack.index( DATA, "preds")$data
    }


    RES = NULL
    RES = spacetime.inla.call( FM=FM, DATA=DATA, SPDE=SPDE, FAMILY="gaussian" )

    # extract summary statistics from a spatial (SPDE) analysis and update the output file
    # inla.summary = spacetime.summary.inla.spde2 = ( RES, SPDE )

    # inla.spde2.matern creates files to disk that are not cleaned up:
    spdetmpfn = SPDE$f$spde2.prefix
    fns = list.files( dirname( spdetmpfn ), all.files=TRUE, full.names=TRUE, recursive=TRUE, include.dirs=TRUE )
    oo = grep( basename(spdetmpfn), fns )
    if(length(oo)>0) file.remove( sort(fns[oo], decreasing=TRUE) )

    rm( SPDE, DATA ); gc()

    # predict upon grid

    if ( method=="direct" ) {
      # direct method ... way too slow to use for production runs
      preds = as.data.frame( plocs )
      preds$xmean = spacetime.invlink( RES$summary.fitted.values[ i_data, "mean"] )
      preds$xsd   = spacetime.invlink( RES$summary.fitted.values[ i_data, "sd"] )
      rm(RES, MESH); gc()
    }

    if (method=="fast") {
      posterior.extract = function(s, rnm) {
        # rnm are the rownames that will contain info about the indices ..
        # optimally the grep search should only be done once but doing so would
        # make it difficult to implement in a simple structure/manner ...
        # the overhead is minimal relative to the speed of modelling and posterior sampling
        i_intercept = grep("intercept", rnm, fixed=TRUE ) # matching the model index "intercept" above .. etc
        i_spatial.field = grep("spatial.field", rnm, fixed=TRUE )
        return(  s$latent[i_intercept,1] + s$latent[ i_spatial.field,1] )
      }

      # note: plocs seems to be treated as token locations and really its range and dims controls output
      pG = inla.mesh.projector( MESH, loc=as.matrix( plocs ) )
      posterior.samples = inla.posterior.sample(n=nsamples, RES)
      rm(RES, MESH); gc()

      rnm = rownames(posterior.samples[[1]]$latent )
      posterior = sapply( posterior.samples, posterior.extract, rnm=rnm )
      posterior = spacetime.invlink( posterior )   # return to original scale
      rm(posterior.samples); gc()

          # robustify the predictions by trimming extreme values .. will have minimal effect upon mean
          # but variance estimates should be useful/more stable as the tails are sometimes quite long
          for (ii in 1:nrow(posterior )) {
            qnt = quantile( posterior[ii,], probs=c(0.025, 0.975), na.rm=TRUE )
            toolow = which( posterior[ii,] < qnt[1] )
            toohigh = which (posterior[ii,] > qnt[2] )
            if (length( toolow) > 0 ) posterior[ii,toolow] = qnt[1]
            if (length( toohigh) > 0 ) posterior[ii,toohigh] = qnt[2]
          }

      # posterior projection is imperfect right now .. not matching the actual requested locations
      preds = data.frame( plon=pG$loc[,1], plat = pG$loc[,2])
      preds$xmean = c( inla.mesh.project( pG, field=apply( posterior, 1, mean, na.rm=TRUE )  ))
      preds$xsd   = c( inla.mesh.project( pG, field=apply( posterior, 1, sd, na.rm=TRUE )  ))
      rm (pG)
    }

    if (0) {
      require(lattice)
      levelplot( log( xmean)  ~ plon+plat, preds, aspect="iso",
                labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      levelplot( log (xsd )  ~ plon+plat, preds, aspect="iso",
                labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }

    return( preds )
  }


