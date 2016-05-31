

spacetime.interpolate.spBayes = function( x, y, z, locsout, n.samples=5000, burn.in=0.5, method="gstat" ) {
  #// interpolate via spBayes ~ random field approx via matern
  #// method determines starting parameter est method to seed spBayes
  require(spBayes)
  library(MBA)
  require( geoR )
  require( stat )
  require( sp )
  require( rgdal )

  if (0) {
    # just for debugging / testing ... and example of access method:
    bioLibrary("bio.utilities", "bio.spacetime")
    require(sp)
    data(meuse)
    x = meuse[, c("x")]
    y = meuse[, c("y")]
    z = meuse$zinc
    n.samples=5000
    burn.in=0.5
    method="gstat"
    locsout = expand.grid(
      seq( min(x), max(x), length.out=500),
      seq( min(y), max(y), length.out=500) )
  }

  xy = as.data.frame( cbind( x,y) )
  z = z / sd( z, na.rm=TRUE)

  pCovars = as.matrix( rep(1, nrow(locsout)))  # in the simplest model, 1 col matrix for the intercept

  stv = spacetime.variogram( xy, z, methods=method )
  rbounds = stv[[method]]$range * c( 0.01, 1.5 )
  phibounds = range( -log(0.05) / rbounds ) ## approximate
  nubounds = c(1e-3, stv[[method]]$kappa * 1.5 )# Finley et al 2007 suggest limiting this to (0,2)
  # Finley, Banerjee Carlin suggest that kappa_geoR ( =nu_spBayes ) > 2 are indistinguishable .. identifiability problems cause slow solutions

  starting = list( phi=median(phibounds), sigma.sq=0.5, tau.sq=0.5, nu=1  ) # generic start
  tuning   = list( phi=starting$phi/10, sigma.sq=starting$sigma.sq/10, tau.sq=starting$tau.sq/10, nu=starting$nu/10 ) # MH variance
  priors   = list(
    beta.flat = TRUE,
    phi.unif  = phibounds,
    sigma.sq.ig = c(5, 0.5), # inverse -gamma (shape, scale):: scale identifies centre; shape higher = more centered .. assuming tau ~ sigma
    tau.sq.ig = c(5, 0.5),  # inverse gamma (shape, scale) :: invGamma( 3,1) -> modal peaking < 1, center near 1, long tailed
    nu.unif = nubounds
  )

  model = spLM( z ~ 1, coords=as.matrix(xy), starting=starting, tuning=tuning, priors=priors, cov.model="matern",
    n.samples=n.samples, verbose=TRUE )

  ##recover beta and spatial random effects
  mm <- spRecover(model, start=burn.in*n.samples )
  mm.pred <- spPredict(mm, pred.covars=pCovars, pred.coords=as.matrix(locsout), start=burn.in*n.samples )
  res = apply(mm.pred[["p.y.predictive.samples"]], 2, mean)



  u = apply(mm$p.theta.recover.samples, 2, mean)
  vrange = geoR::practicalRange("matern", phi=1/u["phi"], kappa=u["nu"]  )

  spb = list( model=model, recover=mm,
    range=vrange, varSpatial=u["sigma.sq"], varObs=u["tau.sq"],  phi=1/u["phi"], kappa=u["nu"] )  # output using geoR nomenclature

  if (plotdata) {
    x11()
    # to plot variogram
    x = seq( 0, vrange* 2, length.out=100 )
    acor = geoR::matern( x, phi=1/u["phi"], kappa=u["nu"] )
    acov = u["tau.sq"] +  u["sigma.sq"] * (1- acor)  ## geoR is 1/2 of gstat and RandomFields gamma's
    plot( acov ~ x , col="orange", type="l", lwd=2, ylim=c(0,max(acov)*1.1) )
    abline( h=u["tau.sq"] + u["sigma.sq"]  )
    abline( h=u["tau.sq"] )
    abline( h=0 )
    abline( v=0 )
    abline( v=vrange )

    round(summary(mm$p.theta.recover.samples)$quantiles,2)
    round(summary(mm$p.beta.recover.samples)$quantiles,2)
    mm.w.summary <- summary(mcmc(t(mm$p.w.recover.samples)))$quantiles[,c(3,1,5)]

    plot(z, mm.w.summary[,1], xlab="Observed w", ylab="Fitted w",
        xlim=range(w), ylim=range(mm.w.summary), main="Spatial random effects")
    arrows(z, mm.w.summary[,1], w, mm.w.summary[,2], length=0.02, angle=90)
    arrows(z, mm.w.summary[,1], w, mm.w.summary[,3], length=0.02, angle=90)
    lines(range(z), range(z))

    obs.surf <-   mba.surf(cbind(xy, z), no.X=500, no.Y=500, extend=T)$xyz.est
    image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
    points(xy)
    contour(obs.surf, add=T)

    x11()
    require(lattice)
    levelplot( res ~ locsout[,1] + locsout[,2], add=T )


  }

  return( res )
}




