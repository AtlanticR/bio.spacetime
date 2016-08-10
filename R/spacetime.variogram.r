
spacetime.variogram = function( xy, z, plotdata=FALSE, edge=c(1/3, 1), methods=c("geoR"), maxdist=NA, nbreaks = 15, functionalform="matern" ) {

  #\\ estimate empirical variograms (actually correlation functions) and then model them using a number of different approaches .. mostly using Matern as basis
  #\\ returns empirical variogram and parameter estimates, and the models themselves
  #\\ expect xy = c(p/lon, p/lat), z= variable
  #\\ ---> removed--> varZ is the total variance which needs to be mulitplied to the curve if you want the "true" semivariance
  #\\ NOTE:: the default parameterization is as in spBayes and gstat which is:

  #\\ matern covariogram (||x||) = sigma^2 * (2^{nu-1} * Gamma(nu) )^{-1} * (phi*||x||)^{nu} * K_{nu}(phi*||x||)
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is the range parameter  
  # -------------------------


  if ( 0 ) {
   # just for debugging / testing ... and example of access method:
   bioLibrary("bio.utilities", "bio.spacetime")
   require(sp)
   data(meuse)
    xy = meuse[, c("x", "y")]
    mz = log( meuse$zinc )
    mm = lm( mz ~ sqrt( meuse$dist ) )
    z = residuals( mm)
   
    plotdata=TRUE
    maxdist =800
    edge=c(1/3, 1)
    nbreaks = 15

    methods=c("gstat", "inla", "geoR" )

    # out = spacetime.variogram( xy, z, methods="spBayes" )
    out = spacetime.variogram( xy, z, methods="spBayes" )
    nd = nrow(out$spBayes$recover$p.theta.samples)
    rr = rep(NA, nd )
    for (i in 1:nd) rr[i] = geoR::practicalRange("matern", phi=1/out$spBayes$recover$p.theta.samples[i,3], kappa=out$spBayes$recover$p.theta.samples[i,4] )
    hist(rr)  # range estimate

    hist( out$spBayes$recover$p.theta.samples[,1] ) #"sigma.sq"
    hist( out$spBayes$recover$p.theta.samples[,2] ) # "tau.sq"
    hist( out$spBayes$recover$p.theta.samples[,3] ) # 1/phi
    hist( out$spBayes$recover$p.theta.samples[,4] ) # nu

    gr = spacetime.variogram( xy, z, methods="geoR" )
    gs = spacetime.variogram( xy, z, methods="gstat" )
    grf = spacetime.variogram( xy, z, methods="RandomFields" )
    gsp = spacetime.variogram( xy, z, methods="spBayes" )
    ginla = spacetime.variogram( xy, z, methods="inla" )

    # tests:

    out = spacetime.variogram( xy, z )
    (out$geoR$range)
    out = spacetime.variogram( xy, z, nbreaks=30 )
    (out$geoR$range)

    out = spacetime.variogram( xy, log(z), nbreaks=30 )
    (out$geoR$range)
    out = spacetime.variogram( xy, log(z) )
    (out$geoR$range)
    require(mgcv)
    og = gam( log(z) ~ s( x) + s(y) + s(x,y), data=xy )
    zr = residuals(og)
    out = spacetime.variogram( xy, zr )  # remove spatial trend results in no variogram, as would be expected
    (out$geoR$range)
    og = gam( log(z) ~ s( elev ) , data=meuse )
    zr = residuals(og)
    out = spacetime.variogram( xy, zr )  # remove spatial trend results in no variogram, as would be expected
    (out$geoR$range)

    require(geoR)
    # plot( out$geoR$vgm )
    # lines( out$geoR$fit, lwd=2, col="slateblue" )
    xRange = c( 0, max(out$geoR$range*2.1 ) )
    yRange = c( 0, max(out$geoR$vgm$v*out$varZ )*1.05 )
    plot ( out$varZ * out$geoR$vgm$v ~ out$geoR$vgm$u, pch=20, xlim=xRange, ylim=yRange, ylab="Semivariance", xlab="Distance" )
      abline( h=0,  col="gray", lwd=2 )
      abline( h= out$varZ *(out$geoR$varSpatial + out$geoR$varObs), lty="dashed", col="slategray"  )
      abline( h= out$varZ * out$geoR$varObs , lty="dashed", col="slategray")
      abline( v=out$geoR$range, lty="dotted", col="slateblue" )
      abline( v=0,  col="gray", lwd=2 )
      x = seq( 0, 2*out$geoR$range, length.out=100 )
      acor = geoR::matern( x, phi=out$geoR$phi, kappa=out$geoR$kappa  )
      acov = out$geoR$varObs +  out$geoR$varSpatial * (1- acor)
      lines( out$varZ * acov ~ x , col="blue", lwd=2 )
  }

  nc_max = 5  # max number of iterations

  out = list()
  out$varZ = var( z, na.rm=TRUE )  # this is the scaling factor for semivariance .. diving by sd, below reduces numerical floating point issues
  out$meanZ = mean(z, na.rm=TRUE)
  out$meanX = mean( xy[,1], na.rm=TRUE )
  out$meanY = mean( xy[,2], na.rm=TRUE )

  z = (z - out$meanZ )/ sqrt( out$varZ ) # (centered and scaled by sd to avoid floating point issues)
  zrange = range( z, na.rm=TRUE )

  names(xy) =  c("plon", "plat" ) # arbitrary
  xr = range( xy$plon, na.rm=TRUE )
  yr = range( xy$plat, na.rm=TRUE )
  drange = min( diff( xr), diff( yr)  )

  # if max dist not given, make a sensible choice
  if ( is.na(maxdist)) {
    maxdist = drange * 0.5  # default
  } else if ( maxdist=="all") {
    maxdist = drange
  }

  xy$plon = ( xy$plon - out$meanX ) / drange
  xy$plat = ( xy$plat - out$meanY ) / drange

  xrange = range( xy$plon, na.rm=TRUE )
  yrange = range( xy$plat, na.rm=TRUE )

  difx = diff( xrange)
  dify = diff( yrange)

  nn = 400
  nxout = trunc(nn * difx / dify)
  nyout = nn
  nzout = 100

  xx = seq( xrange[1], xrange[2], length.out=nxout )
  yy = seq( yrange[1], yrange[2], length.out=nyout )
  zz = seq( zrange[1], zrange[2], length.out=nzout )
  preds = expand.grid( plon=xx, plat=yy )


  # ------------------------

  if ("gstat" %in% methods){
    #\\ covariogram (||x||) = tau^2 * (2^{nu-1} * Gamma(nu) )^{-1} * (phi*||x||)^{nu} * K_{nu}(phi*||x||)
    #\\ gstat::kappa == spBayes::nu
    #\\ gstat::range == spBayes::phi {aka, "scale parameter"}

    require(gstat)
    require(sp)
    vrange = maxdist/2 + 1
    maxdist = maxdist/2 ## back it up a bit to enter smoothly into the loop
    nc = 0
    while ( (maxdist-vrange) < vrange/2  ) {
      nc = nc  + 1
      maxdist = maxdist * 1.5
      vEm = try( variogram( z~1, locations=~plon+plat, data=xy, cutoff=maxdist, width=maxdist/nbreaks, cressie=TRUE ) ) # empirical variogram
      if  ("try-error" %in% vEm) return(NULL)
      vMod0 = vgm(psill=0.75, model="Mat", range=maxdist, nugget=0.25, kappa=1 ) # starting model parameters
      #vMod0 = vgm("Mat")
      vFitgs =  try( fit.variogram( vEm, vMod0, fit.kappa =TRUE, fit.sills=TRUE, fit.ranges=TRUE ) ) ## gstat's kappa is the Bessel function's "nu" smoothness parameter
      if  ("try-error" %in% vFitgs) return(NULL)
      out$gstat = list( fit=vFitgs, vgm=vEm, range=NA, nu=vFitgs$kappa[2], phi=vFitgs$range[2],
          varSpatial=vFitgs$psill[2], varObs=vFitgs$psill[1]  )  # gstat::"range" == range parameter == phi
      out$gstat$range =  geoR::practicalRange("matern", phi=out$gstat$phi, kappa=out$gstat$nu  )
      vrange = out$gstat$range
      if (nc > nc_max ) break()
    }

    if (plotdata) {
      x11()
      plot(vEm, model=vFitgs, add=T)
      x11()
      plot( gamma ~ dist, data=out$gstat$vgm, ylim=c(0,max(out$gstat$vgm$gamma)*1.1), col="blue", pch=20 )
      abline( h=out$gstat$varSpatial + out$gstat$varObs  )
      abline( h=out$gstat$varObs )
      abline( v=out$gstat$range )
      x = seq( 0, maxdist, length.out=100 )
      acor = geoR::matern( x, phi=out$gstat$phi, kappa=out$gstat$nu  )
      acov = out$gstat$varObs + out$gstat$varSpatial * (1- acor)
      lines( acov~x , col="red" )

      if (0) {
        g <- gstat(id = "elev", formula = z~1, locations = ~plon+plat, data = xy )
        g = gstat(g, id="elev", model=vFitgs)
        gpredres <- predict( g, preds )
        x11()
        lp = levelplot( elev.pred ~ plon+plat, gpredres, aspect = "iso", at=zz, col.regions=color.code( "seis", zz),
          contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE)  )
        plot(lp)
      }
    }
  }


  # -------------------------

  if ("geoR" %in% methods) {
    # weighted least squares
    #  he Matern model (correlation function, rho) is defined as:
    #  rho(u;phi,kappa) =(2^(kappa-1) Gamma(kappa))^(-1) (u/phi)^kappa K_kappa(u/phi)
     # where phi and kappa are parameters and K_kappa(...) denotes the
     # modified Bessel function of the third kind of order kappa.  The
     # family is valid for phi > 0 and kappa > 0.
     #\\ default covariogram (||x||) = tau^2 * (2^{nu-1} * Gamma(nu) )^{-1} * (phi*||x||)^{nu} * K_{nu}(phi*||x||)
      #\\ geoR:: rho(h) = (1/(2^(kappa-1) * Gamma(kappa))) * ((h/phi)^kappa) * K_{kappa}(h/phi)

    require( geoR )
    vrange = maxdist/2 + 1
    maxdist = maxdist/2 ## back it up a bit to enter smoothly into the loop
    nc = 0
    while ( (maxdist-vrange) < vrange/2  ) {
      nc = nc + 1
      maxdist = maxdist * 1.25
      vEm = try( variog( coords=xy, data=z, uvec=nbreaks, max.dist=maxdist ) )
      if  ("try-error" %in% vEm) return(NULL)
      vMod = try( variofit( vEm, nugget=0.5, kappa=1, cov.model="matern", ini.cov.pars=c(0.5, maxdist/4) ,
        fix.kappa=FALSE, fix.nugget=FALSE, max.dist=maxdist, weights="cressie" ) )
        # kappa is the smoothness parameter , also called "nu" by others incl. RF
      if  ("try-error" %in% vMod) return(NULL)
      # maximum likelihood method does not work well with Matern
      ML = FALSE
      if (ML) {
        vMod = likfit( coords=xy, data=z, cov.model="matern", ini.cov.pars=vMod$cov.pars,
        fix.kappa=FALSE, fix.nugget=FALSE, kappa=vMod$kappa, nugget=vMod$nugget, lik.method = "REML" )
      }
     vrange = vMod$practicalRange
     if (nc > nc_max ) break()
    }
    out$geoR = list( fit=vMod, vgm=vEm, model=vMod, range=vMod$practicalRange,
              varSpatial= vMod$cov.pars[1], varObs=vMod$nugget, nu=vMod$kappa,  phi=vMod$cov.pars[2] )

    if (plotdata) {
      x11()
      plot(vEm)
      lines(vMod)
      x11()
      plot( out$geoR$vgm )
      points( out$geoR$vgm$v ~ out$geoR$vgm$u, pch=20 )
      abline( h=out$geoR$varSpatial + out$geoR$varObs  )
      abline( h=out$geoR$varObs )
      abline( v=out$geoR$range )
      x = seq( 0, max(out$geoR$vgm$u), length.out=100 )
      acor = geoR::matern( x, phi=out$geoR$phi, kappa=out$geoR$nu  )
      acov = out$geoR$varObs +  out$geoR$varSpatial * (1-acor)  ## geoR is 1/2 of gstat and RandomFields gamma's
      lines( acov ~ x , col="orange" )
    }
  }


  # -------------------------


  if ("RandomFields" %in% methods) {
    require( RandomFields ) ## max likilihood
    rownames( xy) = 1:nrow(xy)  # seems to require rownames ...
    rfdata <- RFspatialPointsDataFrame(
      coords = xy,
        data = z,
        RFparams=list(vdim=1, n=1)
    )
   # uses all data
   # where nu > 0 and K_nu is the modified Bessel function of second kind and distance r >= 0 between two pointsd
   # The Matern covariance model is given by: C(h) = v * phi(A*h/s).
   #  Cov(r) = 2^{1- nu} Gamma(nu)^{-1} (sqrt{2nu} r)^nu K_nu(sqrt{2nu} r)
   # phi = sqrt{2nu}
   
   # RFoptions(
   #   allowdistanceZero=TRUE,
    #  modus_operandi="precise", #‘"careless"’,‘"sloppy"’, ‘"easygoing"’, ‘"normal"’, ‘"precise"’,        ‘"pedantic"’, ‘"neurotic"’
   #   bin_dist_factor=maxdist/2,
      #bins=nbreaks,
      #critical=TRUE, 
   #   approx_zero=0.05, #  Value below which a correlation is considered to be essentially zero.
   #   spConform=TRUE # FALSE is faster
    #)

    model = ~ RMmatern( nu=NA, var=NA, scale=NA) + RMnugget(var=NA)
    
    o = RFfit(model, data=rfdata )
    oo=summary(o)

    out$RandomFields = list ( fit=o, vgm=o[2], model=oo, range=NA,
              varSpatial=oo$param["value", "matern.var"],
              varObs=oo$param["value", "nugget.var"],
              phi=oo$param["value", "matern.s"],  # sqrt(2*RF::".s = scale") == geoR::phi -- need to confirm this :: confirmed JC Jul 2016
              nu=oo$param["value", "matern.nu"], # RF::nu == geoR:: kappa (bessel smoothness param)
              error=NA )

    out$RandomFields$range = geoR::practicalRange("matern", phi=out$RandomFields$phi, kappa=out$RandomFields$nu  )

    if (plotdata) {
      x11()
      plot(  out$RandomFields$vgm@emp.vario ~ out$RandomFields$vgm@centers, pch=20, ylim=c(0,var(z)*1.25) )
      abline( h=out$RandomFields$varSpatial + out$RandomFields$varObs  )
      abline( h=out$RandomFields$varObs )
      abline( v=out$RandomFields$range )

      x = seq( 0, max(out$RandomFields$vgm@centers), length.out=100 )
      acor = geoR::matern( x, phi=out$RandomFields$phi, kappa=out$RandomFields$nu  )
      acov = out$RandomFields$varObs  +  out$RandomFields$varSpatial*(1- acor)
      lines( acov~x , col="red" )

      # compare with:
      x11()
      plot(o)
    }

  }

  # -------------------------

  if ("spBayes" %in% methods) {
    # note spBayes::phi = 1/ gstat::phi  
    require(spBayes)
    library(MBA)
    require( geoR )
    geoR = spacetime.variogram( xy, z, methods="geoR" )
    rbounds = c( median( diff(  geoR$geoR$vgm$u) )/2, geoR$geoR$range *1.5 )
    phibounds = range( -log(0.05) / rbounds ) ## approximate
    nubounds = c(1e-3, geoR$geoR$nu * 1.5 )# Finley et al 2007 suggest limiting this to (0,2)
    # Finley, Banerjee Carlin suggest that kappa_geoR ( =nu_spBayes ) > 2 are indistinguishable .. identifiability problems cause slow solutions
    n.samples = 5000
    starting = list( phi=median(phibounds), sigma.sq=0.51, tau.sq=0.51, nu=1.1  ) # generic start
    #starting = list( phi=1/2, sigma.sq=res.geoR$geoR$varSpatial, tau.sq=res.geoR$geoR$varObs, nu=30  ) # generic start
    tuning   = list( phi=starting$phi/12, sigma.sq=starting$sigma.sq/12, tau.sq=starting$tau.sq/12, nu=starting$nu/12 ) # MH variance to get acceptance rante bet 30-40%
    priors   = list(
      beta.flat = TRUE,
      phi.unif  = phibounds,
      sigma.sq.ig = c(5, 0.5),  # inverse -gamma (shape, scale):: scale identifies centre; shape higher = more centered .. assuming tau ~ sigma
      tau.sq.ig = c(5, 0.5),  # inverse gamma (shape, scale) :: invGamma( 3,1) -> modal peaking < 1, center near 1, long tailed
      nu.unif = nubounds
    )

    model = spLM( z ~ 1, coords=as.matrix(xy), starting=starting, tuning=tuning, priors=priors, cov.model="matern",
      n.samples=n.samples, verbose=TRUE )

    burn.in <- 0.5*n.samples

    ##recover beta and spatial random effects
    m.1 <- spRecover(model, start=burn.in )

    u = apply(m.1$p.theta.recover.samples, 2, mean)
    vrange = geoR::practicalRange("matern", phi=1/u["phi"], kappa=u["nu"]  )

    out$spBayes = list( model=model, recover=m.1,
      range=vrange, varSpatial=u["sigma.sq"], varObs=u["tau.sq"],  phi=1/u["phi"], nu=u["nu"] )  # output using geoR nomenclature

    if (plotdata) {
      x11()
      x = seq( 0, vrange* 2, length.out=100 )
      acor = geoR::matern( x, phi=1/u["phi"], kappa=u["nu"] )
      acov = u["tau.sq"] +  u["sigma.sq"] * (1- acor)  ## geoR is 1/2 of gstat and RandomFields gamma's
      plot( acov ~ x , col="orange", type="l", lwd=2, ylim=c(0,max(acov)*1.1) )
      abline( h=u["tau.sq"] + u["sigma.sq"]  )
      abline( h=u["tau.sq"] )
      abline( h=0 )
      abline( v=0 )
      abline( v=vrange )

      if(0) {
        round(summary(m.1$p.theta.recover.samples)$quantiles,2)
        round(summary(m.1$p.beta.recover.samples)$quantiles,2)
        m.1.w.summary <- summary(mcmc(t(m.1$p.w.recover.samples)))$quantiles[,c(3,1,5)]

        plot(z, m.1.w.summary[,1], xlab="Observed w", ylab="Fitted w",
            xlim=range(z), ylim=range(m.1.w.summary), main="Spatial random effects")
        arrows(z, m.1.w.summary[,1], z, m.1.w.summary[,2], length=0.02, angle=90)
        arrows(z, m.1.w.summary[,1], z, m.1.w.summary[,3], length=0.02, angle=90)
        lines(range(z), range(z))

        par(mfrow=c(1,2))
        obs.surf <-   mba.surf(cbind(xy, z), no.X=100, no.Y=100, extend=T)$xyz.est
        image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
        points(xy)
        contour(obs.surf, add=T)
      }
    }

  }


  # -------------------------


  if ("inla" %in% methods){
    require(INLA)
    require(lattice)

    inla.setOption(scale.model.default = TRUE)  # better numerical performance of IGMRF models and less dependnence upon hyperpriors

    locs0  = as.matrix( xy )
    xy$b0 = 1  # intercept for inla

    vRange = maxdist/10

    M0.domain = inla.nonconvex.hull( locs0 )
    MESH = inla.mesh.2d (
      loc=locs0, # locations of data points
      boundary = M0.domain,
      max.edge = edge * vRange
    )

#    kappa0 = sqrt(8) / vRange
#    tau0 = 1/ ( sqrt(4*pi) * kappa0 * vPsill )
    alpha = 2 # -> alpha-1 == nu (inla fixes it at 1,2,or 3)
    SPDE = inla.spde2.matern( MESH, alpha=alpha )
    spatial.field <- inla.spde.make.index('spatial.field', n.spde=SPDE$n.spde )

    # projection matrix A to translate from mesh nodes to data nodes
    A = inla.spde.make.A( mesh=MESH, loc=locs0 )

    # data stack for occurence (PA)
    Z = inla.stack(
        tag="data",
        data=list( z=z ) ,
        A=list(A, 1 ),
        effects=list( spatial.field=spatial.field, xy )  # b0 is the intercept
    )

    RES <- inla(  z ~ 0 + b0+ f( spatial.field, model=SPDE ), family="gaussian",
        data=inla.stack.data(Z),
        control.compute=list(dic=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.fixed = list(expand.factor.strategy='inla') ,
        control.predictor=list(A=inla.stack.A(Z), compute=TRUE, link=1 ) ,
        control.inla = list( h=1e-4, tolerance=1e-10),
        # control.inla=list(strategy="laplace", npoints=21, stencil=7 ) ,
        verbose = FALSE
    )

    oo = inla.spde2.result(RES, "spatial.field", SPDE, do.transf=TRUE)

    inames = c( "mode", "mean", "sd", "quant0.025", "quant0.25", "quant0.5",  "quant0.75", "quant0.975", "low", "high" )

    # Range parameter .. ie, sqrt(8)/exp(oo$summary.log.kappa$mean)
    im = oo$marginals.range.nominal[[1]]
    iRange = c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im )) )

    # "Spatial variance/error ('partial sill variance')"
    im = oo$marginals.variance.nominal[[1]]
    iVar =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im )) )

    # kappa.inla  == 1/phi.geoR
    im = oo$marginals.kappa[[1]]
    iKappa =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im ) ) )

    # tau
    im = oo$marginals.tau[[1]]
    iTau =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im ) ) )

    ## Non-spatial ("observation") error ('nugget variance')
    iprec = grep ( "Precision.*observ.*", names(RES$marginals.hyperpar), ignore.case=TRUE )
    im = inla.tmarginal( function(x) {1/x}, RES$marginals.hyperpar[[ iprec ]] )
    iNugget =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im ) ) )

    inla.summary = as.matrix( rbind( iKappa, iTau, iRange, iVar, iNugget ) )
    rownames( inla.summary) = c( "kappa", "tau", "range", "spatial error", "observation error" )
    colnames( inla.summary) = inames

    out$inla = list(mesh=MESH, res=RES, range.inla90=inla.summary[["range","mean"]],
      varSpatial=inla.summary[["spatial error","mean"]], varObs=inla.summary[["observation error","mean"]],
      phi = 1/inla.summary[["kappa","mean"]] , nu=alpha-1, error=NA )

    # kappa{geoR} = lambda{INLA} == alpha-1 {INLA} and alpha=2 by default in INLA
    out$inla$range = geoR::practicalRange("matern", phi=out$inla$phi, kappa=out$inla$nu  )


    if (plotdata) {
      require( geoR )
      x = seq( 0,  out$inla$range * 1.5, length.out=100 )
      svar =  out$inla$varObs + out$inla$varSpatial * (1-geoR::matern( x, phi=out$inla$phi, kappa=out$inla$nu  ))
      plot( svar~x, type="l" )
      abline( h=out$inla$varObs + out$inla$varSpatial )
      abline( h=out$inla$varObs )
      abline( v=out$inla$range.inla90  )
      abline( v=out$inla$range, col="red"  )
    }

  }


  # -------------------------

  if ("BayesX" %in% methods){
    library("R2BayesX")
    # fixes nu=1.5
    # phi = max(distance) / const, such that Corr(distance=const) = 0.001; i.e. range at distance where covar ~0.999
    # not sure how to recover the correct phi/range from this ...

    fm1 <- bayesx( z ~ sx(plon, plat, bs="kr" ), family="gaussian", method="REML", data =xy )
    #out$BayesX = list( fit=vMod, vgm=vEm, model=vMod, range=vMod$practicalRange,
    #          varSpatial= vMod$cov.pars[1], varObs=vMod$nugget, nu=1.5,  phi=vMod$cov.pars[2] )
    summary( fm1)
        # Call:
        # bayesx(formula = z ~ sx(plon, plat, bs = "kr"), data = xy, family = "gaussian", 
        #     method = "REML")
         
        # Fixed effects estimation results:

        # Parametric coefficients:
        #             Estimate Std. Error t value  Pr(>|t|)    
        # (Intercept)   8.1534     0.4787  17.033 < 2.2e-16 ***
        # ---
        # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

        # Smooth terms:
        #               Variance Smooth Par.     df Stopped
        # sx(plon,plat)   1.5145      0.4342 8.4592       0
         
        # Scale estimate: 0.6577 
         
        # N = 155  df = 9.45917  AIC = 99.5061  BIC = 128.294  
        # GCV = 0.700414  logLik = -40.2939  method = REML  family = gaussian  
      return(fm1)
  }


  # -------------------------

  if ("jags" %in% methods){
    require(rjags)
    require(jagsUI)
    # assume nu = 1 (due to identifiability issues)

    print( "Slow ... 10 min for meuse test data")

    jagsmodel = paste0("
    model{
      for(i in 1:N){
        y[i] ~ dnorm(mu[i], prec)
        mu[i] = beta0 + errorSpatial[i]
        muSpatial[i] = 0
      }
      prec = 1.0/ (tausq + sigmasq )
      invCOVAR = inverse(COVAR)
      errorSpatial ~ dmnorm( muSpatial, invCOVAR)
      for(i in 1:N) {
        COVAR[i,i] = sigmasq
        for(j in 1:(i-1)) {
          COVAR[i,j] = sigmasq * exp(-( DIST[i,j]/phi))
          COVAR[j,i] = COVAR[i,j]
        } 
      }
      tausq = 1/tausq.inv
      tausq.inv ~ dgamma(0.1,0.1)
      sigmasq = 1/sigmasq.inv
      sigmasq.inv ~ dgamma(2,1)
      phi ~ dgamma(1,0.1)
      beta0 ~ dnorm(0,0.0001)
    } ")

  fn = tempfile()
  cat( jagsmodel, file=fn )
  distances = as.matrix(dist( xy, diag=TRUE, upper=TRUE))
  Data = list( N=length(z), DIST=distances, y=z )
  fit = jagsUI::jags(data=Data, 
       parameters.to.save=c("phi", "sigmasq", "tausq"),
       model.file=fn,
       n.iter=2000,
       n.chains=3,
       n.burnin=200,
       n.thin=10,
       parallel=TRUE,
       DIC=FALSE)
   summary(fit)
   plot(fit)
   gelman.plot(fit$samples)
  # geweke.plot(fit$samples)
  #update(fit, n.iter=2000, n.thin=20 )
    acf( fit$sims.list$phi)
    acf( fit$sims.list$sigmasq)
    acf( fit$sims.list$tausq)

   #  apply( fit$phi, 1, quantile, probs=c(0.025, 0.5, 0.975) )
   return(fit$summary)
  }


  # -------------------------

  if ("TMB" %in% methods){
   
  }

  # -------------------------

  if ("LaplacesDemon" %in% methods){
    
    stop("This is too slow to use. It is just to document the approach. dmvn is the culprit .. try sparse matrix methods ...")
    # uses the gstat::matern parameterization
    # nu is "scale" (hard to identify due to exp(phi*x)^nu == exp(phi*nu*x)) .. maybe just do phi*nu ?

    require(LaplacesDemonCpp)
    Data = list(
      eps = 1e-9,
      N = length(z),  # required for LaplacesDemon
      DIST=as.matrix(dist( xy, diag=TRUE, upper=TRUE)), # distance matrix between knots
      y=z  
    )
    Data$mon.names = c( "LP", paste0("mu",1:Data$N ) )
    Data$parm.names = as.parm.names(list(tau=0, sigma=0, phi=0, nu=0, muSpatial=rep(0,Data$N) ))
    Data$pos = list(
      tau = grep("tau", Data$parm.names),
      sigma = grep("sigma", Data$parm.names),
      phi = grep("phi", Data$parm.names),
      nu = grep("nu", Data$parm.names),
      muSpatial = grep("muSpatial", Data$parm.names)
    )
    Data$PGF = function(Data) {
      tau = rhalfcauchy( 1, 25 )
      sigma = rhalfcauchy( 1, 25 )
      phi = rhalfcauchy( 1, 25 )
      nu = 1
      muSpatial = mvnfast::rmvn(1, rep(0,Data$N), sigma*sigma*exp(-phi*Data$DIST )^nu )
      return(c(muSpatial, tau, sigma, phi, nu))
    }
    Data$PGF  = compiler::cmpfun(Data$PGF)
    Data$Model = function(parm, Data){
      nu = parm[Data$pos$nu] = 1  # in case GIV resets to some other value
      tau = parm[Data$pos$tau] = LaplacesDemonCpp::interval(parm[Data$pos$tau], Data$eps, Inf)
      sigma = parm[Data$pos$sigma] = LaplacesDemonCpp::interval(parm[Data$pos$sigma], Data$eps, Inf)
      phi = parm[Data$pos$phi] = LaplacesDemonCpp::interval(parm[Data$pos$phi], Data$eps, Inf)
      muSpatial = parm[Data$pos$muSpatial]
      covSpatial = sigma*sigma * exp(-phi * Data$DIST)^nu   ## spatial correlation
      muSpatial.prior =  mvnfast::dmvn( muSpatial, rep(0, Data$N), sigma=covSpatial, log=TRUE )
      tau.prior = dgamma(tau, 1, 0.001, log=TRUE)
      sigma.prior = dgamma(sigma, 1, 0.001, log=TRUE)
      phi.prior = dgamma(phi, 1, 0.001, log=TRUE)
      #nu.prior = dgamma(nu, 1, 0.001 log=TRUE)
      nu.prior = 0
      nugget = rnorm(Data$N, 0, tau) 
      yhat =  nugget + muSpatial # local iid error + spatial error
      LL = sum(dnorm(Data$y, yhat, sqrt( tau*tau+sigma*sigma), log=TRUE)) ## Log Likelihood
      LP = LL + muSpatial.prior + sigma.prior + phi.prior + nu.prior ### Log-Posterior
      Modelout = list(LP=LP, Dev=-2*LL, Monitor=c(LP, mu), yhat=mu, parm=parm)
      return(Modelout)
    }

    Data$Model.ML  = compiler::cmpfun( function(...) (Data$Model(...)$Dev / 2) )  # i.e. - log likelihood
    Data$Model.PML = compiler::cmpfun( function(...) (- Data$Model(...)$LP) ) #i.e., - log posterior 
    Data$Model = compiler::cmpfun(Data$Model) #  byte-compiling for more speed .. use RCPP if you want more speed

    print (Data$Model( parm=Data$PGF(Data), Data ) ) # test to see if return values are sensible

    parm0 = Data$PGF(Data)

# maximum likelihood solution
f.ml = optim( par=Data$PGF(Data), fn=Data$Model.ML, Data=Data, control=list(maxit=5000, trace=1), method="BFGS"  )
names(f.ml$par ) = Data$parm.names

# penalized maximum likelihood .. better but still a little unstable depending on algorithm
f.pml = optim( par=Data$PGF(Data), fn=Data$Model.PML, Data=Data,  control=list(maxit=5000, trace=1), method="BFGS" , hessian=FALSE )
names(f.pml$par ) = Data$parm.names
#print(sqrt( diag( solve(f.pml$hessian) )) ) # assymptotic standard errors


    # burn-in
    f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=parm0, Iterations=100, Status=10, Thinning=10 )
    plot(f, Data=Data)

    # sampling
    f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=1000, Status=10, Thinning=1, Covar=f$Covar )
   
    f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=1000, Status=10, Thinning=1, Covar=f$Covar, Method="NUTS" )
   
    # first mcmc .. faster start-up and convergence to global equil 
    f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=Data$PGF(Data), Iterations=500, Status=1, Thinning=1 )
    
    # now use the currest parm est and try a laplace:
    f = LaplaceApproximation(Data$Model, Data=Data, parm=Data$PGF(Data), Iterations=500, Method="Roptim", CPUs=6, method="BFGS" ) # refine it

    f = LaplaceApproximation(Data$Model, Data=Data, parm=as.initial.values(f), Iterations=50, Method="Roptim", CPUs=6, method="BFGS" ) # refine it
    
    # f = VariationalBayes(Data$Model, Data=Data, parm=as.initial.values(f), Iterations=500, Samples=20, CPUs=5, Covar=f$Covar ) # refine it again
    # f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=1000, Status=10, Thinning=1, Covar=f$Covar )
    # f = VariationalBayes(Data$Model, Data=Data, parm=as.initial.values(f), Iterations=100, Samples=10, CPUs=5, Covar=f$Covar )
    # f = IterativeQuadrature(Data$Model, Data=Data, parm=as.initial.values(f), Iterations=10, Algorithm="AGH",
    #  Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL), Covar=f$Covar )
        
   f = LaplaceApproximation(Data$Model, Data=Data, parm=as.initial.values(f), Method="TR", Iterations=1000 , CPUs=6 ) 
    f = LaplaceApproximation(Data$Model, Data=Data, parm=as.initial.values(f), Method="BFGS", Iterations=1000, CPUs=6  ) 
     f = LaplaceApproximation(Data$Model, Data=Data, parm=as.initial.values(f), Method="SPG", Iterations=1000, CPUs=6  ) 

    if (plotdata) {
      Consort(f)
      plot(f)
       m = f$Summary2[grep( "(\\<mu)([0123456789]{1,3})", rownames( f$Summary2 ) ),]
      # m = f$Summary2[grep( "muSpatial", rownames( f$Summary2 ) ),]
      plot( Data$y~m[, "Mean"], pch="." )
    }
  }

  out$LaplacesDemon = list( fit=f, vgm=NA, model=NA, range=vMod$practicalRange,
            varSpatial= vMod$cov.pars[1], varObs=vMod$nugget, nu=1.5,  phi=vMod$cov.pars[2] )
  out$gstat$range =  geoR::practicalRange("matern", phi=out$gstat$phi, kappa=out$gstat$nu  )

  return(out)
}


