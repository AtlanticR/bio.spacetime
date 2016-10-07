
  # example code to use the spacetime interpolation module
  # inputs: lon, lat, time, var as an array
  # internally converts to planar coords and then grids the data

  p = list()
	p$libs = bioLibrary( "bio.utilities", "bio.spacetime", "bio.temperature" )
  p$libs = c( p$libs, RLibrary( "lubridate", "parallel", "sp", "mgcv", "lattice", "gstat", "INLA" )   )
  p$clusters = rep("localhost", detectCores() )

  # p$project.name = "satellite.color"
  p$project.name = "temperature"

  p$project.outdir.root = project.datadirectory( p$project.name, "analysis" )

  p = spacetime_parameters( p, "SSE" )  # data  domain
  p$crs = lookup.projection.params( p$internal.projection )

  p$vars.required = c("timestamp", "longitude", "latitude" )
  p$vars.focal = "temperature"
  p$vars.covar = c( "depth", "salinity", "sigmat", "oxyml" )
  p$vars.keep = c( p$vars.required, p$vars.focal, p$vars.covar )

  # 1. discretise to space, time  blocks or not ..
  p$ti.offset = 0  # to do a phase shift of time
  p$ti.scale = 1/4  # fraction of year to discretize: (i.e., 1 / no time increments in a year )
  p$sp.scale = 10   # spatal scale to discretize: km
  p$edge = c(1/3, 1)

  p$tyears = c(1950:2014)

  V = NULL
  for( yr in p$tyears ) {
    Z = hydro.db( DS="bottom.annual", p=p, yr=yr )
    Z$timestamp = as.POSIXct( Z$date)
    Z = Z[, p$vars.keep]
    V = rbind( V, Z )
  }
  rm (Z)

  V = V[ which( is.finite( V$lon + V$lat ) ) , ]

# fix spatial range and variable names
  uu = which( names(V) == "longitude"); names(V)[uu] = "lon"
  uu = which( names(V) == "latitude");  names(V)[uu] = "lat"

  V = lonlat2planar( V, proj.type=p$crs ) # obtain planar coords

  # limit to area of interest:
  aoi = which( V$lat >= p$corners$lat[1] & V$lat <= p$corners$lat[2] &
               V$lon >= p$corners$lon[1] & V$lon <= p$corners$lon[2] )# basic rectangle around area
  # aoi = which( V$lat >= 42 & V$lat <= 45 & V$lon >= -62 & V$lon <= -58)# basic rectangle around area
  # aoi = which( V$lat >= 38 & V$lat <= 48 & V$lon >= -71 & V$lon <= -48)# basic rectangle around area
  # aoi = filter.region.polygon( V[, c("lon", "lat")], region="maritimes.region")  # to be defined
  V = V[aoi,]

  # V$lon = V$lat = NULL

  # 2. interpolate in time with small spatial extent  .. constraint length autocor length (and n samp)

  # 3. interpolate in space with small time extent .. constraint temporal autocor lenth (or no. samp)
  # repeat until convergence

  # discretize time
  V$ts = discretize.time( V$timestamp, ti.scale=p$ti.scale )

  # V$plon = grid.internal( V$plon, p$plons )
  # V$plat = grid.internal( V$plat, p$plats )


  nzout = 100
  zrange = range( V[, p$vars.focal], na.rm=TRUE )
  zz = seq( zrange[1], zrange[2], length.out=nzout )

  preds = expand.grid( plon=p$plons, plat=p$plats )

  sset = which(is.finite(V$temperature) & lubridate::year(V$timestamp) %in% c(1990:1995) )
  W = V[sset, ]

  # first pass -- use gstat to obtain faster estimates of variogram parameters to speed up inla
  vEm = variogram( temperature ~1, locations=~plon+plat, data=W ) # empirical variogram

  drange = sqrt(diff(range(p$plons))^2 + diff(range(p$plats))^2) # corner to corner distance
  vMod0 = vgm(psill=0.5, model="Mat", range=drange/10, nugget=0.5, kappa=2 ) # starting model parameters
  vFitgs =  fit.variogram( vEm, vMod0 ) ## gstat's kappa is the Bessel function's "nu" smoothness parameter

  plot( vEm, model=vFitgs, main="gstat Variogram" )

  vRange = vFitgs$range[2] # 95% of total variance
  vTot   = vFitgs$psill[1] + vFitgs$psill[2]
  vPsill = vFitgs$psill[2]
  vNugget = vFitgs$psill[1]

  # basic kriging
  g = gstat(id = "temp", formula = temperature ~ 1, locations = ~plon+plat, data = W )
    g = gstat(g, id="temp", model=vFitgs)
    gpredres <- predict( g, preds )
    x11()
    lp = levelplot( gpredres ~ plon+plat, preds, aspect = "iso", at=zz, col.regions=color.code( "seis", zz),
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE)  )
    plot(lp)


  # GAM
  g = gam( temperature ~ s(plon,plat), data = W )
  gpredres <- predict( g, preds )
    lp = levelplot( gpredres ~ plon+plat, preds, aspect = "iso", at=zz, col.regions=color.code( "seis", zz),
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE)  )
    plot(lp)



  # GAM
  g1 = gam( temperature ~ s(plon,plat) + log(depth) , data = W )
  gpredres <- predict( g1, preds )
    lp = levelplot( gpredres ~ plon+plat, preds, aspect = "iso", at=zz, col.regions=color.code( "seis", zz),
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE)  )
    plot(lp)


  # INLA
  inla.setOption(scale.model.default = TRUE)  # better numerical performance of IGMRF models and less dependnence upon hyperpriors

  locs0  = as.matrix(W[, c("plon", "plat")] )
  z = W[,p$var.focal]  # must live outside W for inla
  # W$z = NULL
  W$b0 = 1  # intercept for inla

  edge = 1/2

  M0.domain = inla.nonconvex.hull( locs0 )
  M0 = inla.mesh.2d (
    loc=locs0, # locations of data points
    boundary = M0.domain,
    max.edge = edge * vRange
  )

  # "prior" specifications
  kappa0 = sqrt(8) / vRange
  tau0 = 1/ ( sqrt(4*pi) * kappa0 * vPsill )

  S0 = inla.spde2.matern( M0, alpha=2,
    B.tau = cbind(log(tau0), -1,1 ),     # parameter basis functions
    B.kappa = cbind( log(kappa0), 0, 1 ), # parameter basis functions
    theta.prior.mean = c(0,  0),    # theta1 controls variance .. vague; theta2 controls range   .. means 0
    theta.prior.prec = c(0.1, 1)  #  precisions are vague for theta1;  for range .. theta2 prec 1 ==> 95% prior prob that range is smaller than domain size
  )

  i <- inla.spde.make.index('i', n.spde=S0$n.spde )

  # projection matrix A to translate from mesh nodes to data nodes
  A = inla.spde.make.A( mesh=M0, loc=locs0 )

  # data stack for occurence (PA)
  Z = inla.stack(
      tag="data",
      data=list( z=z ) ,
      A=list(A, 1 ),
      effects=list( i=i, xyz )  # b0 is the intercept
  )

  R <- inla(  z ~ 0 + b0+ f( i, model=S0, diagonal=1e-2),
      data=inla.stack.data(Z),
      control.compute=list(dic=TRUE),
      control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
      control.predictor=list(A=inla.stack.A(Z), compute=TRUE) ,
      # control.inla=list(strategy="laplace", npoints=21, stencil=7 ) ,
      verbose = FALSE
  )

  # field parameters on user scale
  oo = inla.spde2.result(R, 'i', S0, do.transf=TRUE)

  extract =  c("mean", "sd", "mode", "0.5quant", "0.025quant", "0.975quant")
  iRange = exp( oo$summary.log.range.nominal[ extract] ) # or iRange=sqrt(8)/exp(oo$summary.log.kappa$mean)
  iVar = exp( oo$summary.log.variance.nominal[extract] ) # spatial variance (~ psill)
  iKappa = exp( oo$summary.log.kappa[extract]  )
  iTau = exp(oo$summary.log.tau[extract ] )

  iNuggetmarginals = inla.tmarginal( function(x) {1/x}, R$marginals.hyperpar[["Precision for the Gaussian observations"]] )
  iNugmode = inla.mmarginal( iNuggetmarginals )
  iNugmean = inla.emarginal( function(x) {x}, iNuggetmarginals )
  iNugSD =  sqrt( inla.emarginal( function(x) {x^2}, iNuggetmarginals ) - iNugmean^2)
  iNugQuants = inla.qmarginal( c(0.5, 0.025, 0.975), iNuggetmarginals )

  iNugget = data.frame(cbind( iNugmean, iNugSD, iNugmode, iNugQuants[1], iNugQuants[2], iNugQuants[3] ) )
  names(iNugget) = extract

  # indices for random field at data locations
  idat <- inla.stack.index( Z, 'data')$data
  # correlation between the the posterior mean and the response by
  cor.predict = cor( z, R$summary.linear.predictor$mean[idat])

  if (plot) {
    x11();
    plot(M0, asp=1 ) # visualise mesh

    x11()
    vn = "b0"
    v=R$marginals.fixed[[vn]]
    v =  v[order(v[,1]),]
    plot( v, type="l", xlab="b0 -- intercept", ylab="Density" )

    x11()
    vn = "Precision for the Gaussian observations"  ## i.e, the "nugget" or observation error
    v = R$marginals.hyperpar[[vn]]
    # v = v[order(v[,1]),]
    plot.default( iNuggetmarginals, xlab="Non-spatial observation error ('nugget variance')", type="l", ylab="Density" )
    abline( v=iNugget$mean, lty="dotted" )

    x11();
    vn = "marginals.variance.nominal"  # spatially stuctured variance .. ~ psill
    plot( oo[[vn]][[1]], type='l', xlab="Spatial error ('partial sill variance')", ylab='Density')
    abline(v=iVar$mean, lty="dotted" )

    x11();
    vn = "marginals.kappa"
    plot(oo[[vn]][[1]], type='l', xlab=expression(kappa), ylab='Density')
    abline(v=iKappa$mean, lty="dotted"  )

    x11();
    vn = "marginals.range.nominal"
    plot(oo[[vn]][[1]], type='l', xlab='range nominal', ylab='Density')
    abline(v=iRange$mean, lty="dotted"  )

    x11()
    require(lattice)

    delta = mean( R$summary.linear.predictor$mean[idat]) -  mean(R$summary.random$i$mean)
    pG = inla.mesh.projector( M0, xlim=xrange, ylim=yrange, dims=c(nxout, nyout) )
    out = inla.mesh.project( pG, R$summary.random$i$mean ) # mean
    outdf = as.data.frame.table( out)
    preds$z = outdf[,3] + delta
    datarange = range( preds$z, na.rm=TRUE )
    dr = seq( datarange[1], datarange[2], length.out=150)
    lp = levelplot( z~plon+plat, preds, aspect="iso", main="Posterior mean", at=dr, col.regions=color.code( "seis", dr) ,
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    print(lp)

    pG = inla.mesh.projector( M0, xlim=xrange, ylim=yrange, dims=c(nxout, nyout) )
    out = inla.mesh.project( pG, R$summary.random$i$sd ) # SD
    outdf = as.data.frame.table( out)
    preds$z = outdf[,3]
    datarange = range( preds$z, na.rm=TRUE )
    dr = seq( datarange[1], datarange[2], length.out=150)
    lp = levelplot( z~plon+plat, preds, aspect="iso", main="Posterior SD", at=dr, col.regions=color.code( "seis", dr) ,
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    print(lp)


  }

  inla.summary = rbind( iKappa, iTau, iRange, iVar, iNugget )
  rownames( inla.summary) = c( "kappa", "tau", "range", "spatial error", "observation error" )
  colnames( inla.summary) = c( "mean", "sd", "mode", "median", "q025", "q975" )

  out = list( vario.empirical=vEm, vario.model=vFitgs, vario.range=vRange, vario.psill=vPsill, vario.nugget=vNugget,
              kappa0=kappa0, tau0=tau0, inla.summary=inla.summary, correl=cor.predict
  )





