

spacetime.interpolate.xy = function( interp.method, data, locsout, 
  trimquants=TRUE, trimprobs=c(0.025, 0.975), 
  nr=NULL, nc=NULL, theta=NULL, xwidth=theta*10, ywidth=theta*10, 
  link=NA ) {
  #\\ reshape after interpolating to fit the output resolution 
  
  # trim quaniles in case of extreme values
  if (trimquants) {
    vq = quantile( data$z, probs=trimprobs, na.rm=TRUE )
    ii = which( data$z < vq[1] )
    if ( length(ii)>0) data$z[ii] = vq[1]
    jj = which( data$z > vq[2] )
    if ( length(jj)>0) data$z[jj] = vq[2]
  }
  out = res = NULL

  if (interp.method == "kernel.density") {
    # default :: create a "surface" and reshape to a grid using (gaussian) kernel-based smooth via FFT
    require(fields)
    isurf = fields::interp.surface( data, loc=locsout  )
    rm (data); gc()
    zout = matrix( isurf, nrow=nr, ncol=nc )
    out = fields::image.smooth( zout, theta=theta, xwidth=xwidth, ywidth=ywidth ) 
    return (out$z)
  }

  # ------


spacetime.regrid = function( Z0, L0, L1, p0, p1, method="fast" ) {
  #\\ regrid/reproject from p0 to p1 
  #\\ rgdal calls this "warping"
  
  if (method=="raster") {
    # here Z0 has to be rasterizable or a list of rasters or a brick
    require (raster)
    for (vn in names(Z0)) {
       Z[[vn]] = projectRaster( 
          from = rasterize( Z0, spatial.parameters.to.raster(p0), field=vn, fun=mean), 
          to   = spatial.parameters.to.raster( p1) )
    }
    return (Z)
  }

  if (method=="fast") {
      # extract coords, convert to new coords, interpolate where required an
      M = matrix( NA, nrow=p0$nplons, ncol=p0$nplats) # matrix respresentation of the data in new coord system
      L0$plon = grid.internal( L0$plon, p0$plons ) # ensure correct resolution
      L0$plat = grid.internal( L0$plat, p0$plats )
      L2M = cbind( ( L0$plon-p0$plons[1])/p0$pres + 1, (L0$plat-p0$plats[1])/p0$pres + 1) # row, col indices in matrix form of the new coordinate system
      L1 = planar2lonlat( L1, proj.type=p1$internal.projection )
      L1 = lonlat2planar( L1, proj.type=p0$internal.projection )  # convert lon, lat to old projection     
      M[L2M] = Z0
      Z = fields::interp.surface( list( x=p0$plons, y=p0$plats, z=M ), loc=L1[,c("plon","plat")] )
      ii = which( is.na( Z ) )
      if ( length( ii) > 0 ) {
        # try again ..
        Z[ii] = fields::interp.surface( list( x=p0$plons, y=p0$plats, z=M ), loc=L1[ ii, c("plon","plat")] )
      }
      ii = which( is.na( Z ) )
      if ( length( ii) > 0 ) {
        Zii =  fields::image.smooth( M, dx=p0$pres, dy=p0$pres, wght=p0$wght )$z  
        Z[ii] = Zii[ii]
      }
      return( Z)
  }

  # ----------
  
  if (interp.method == "inla.spde" ) {
    require(INLA)
    FM = formula( "Y ~ -1 + intercept + f( spatial.field, model=SPDE )" )
    Y = data$z
    locs= cbind(data$x, data$y)
    rm (data); gc()
    method="fast"  # "direct" is slower and more accurate
    nsamples=5000 
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
    varY = as.character( FM[2] )
    obs_index = inla.spde.make.index(name="spatial.field", SPDE$n.spde)
    obs_eff = list()
    obs_eff[["spde"]] = c( obs_index, list(intercept=1) )
    obs_A = list( inla.spde.make.A( mesh=MESH, loc=locs[,] ) ) # no effects
    obs_ydata = list()
    obs_ydata[[ varY ]] = spacetime.link ( Y )
    DATA = inla.stack( tag="obs", data=obs_ydata, A=obs_A, effects=obs_eff, remove.unused=FALSE )
    if ( method=="direct") {
      # direct method
      preds_index = inla.spde.make.index( name="spatial.field", SPDE$n.spde)
      preds_eff = list()
      preds_eff[["spde"]] = c( list( intercept=rep(1,MESH$n )),
           inla.spde.make.index(name="spatial.field", n.spde=SPDE$n.spde) )
      ydata = list()
      ydata[[ varY ]] = NA
      Apreds = inla.spde.make.A(MESH, loc=as.matrix( locsout ) )
      PREDS = inla.stack( tag="preds", data=list( Y=NA), A=list(Apreds),
        effects=list( c( list(intercept=rep(1, MESH$n)), inla.spde.make.index( name="spatial.field", MESH$n))) )
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
      preds = as.data.frame( locsout )
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
      # note: locsout seems to be treated as token locations and really its range and dims controls output
      pG = inla.mesh.projector( MESH, loc=as.matrix( locsout ) )
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
    return(preds$xmean) 
  }

}
   
