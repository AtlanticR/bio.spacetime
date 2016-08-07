
spacetime.variogram = function( xy, z, plotdata=FALSE, edge=c(1/3, 1), methods=c("geoR"), maxdist=NA, nbreaks = 15 ) {

  #\\ estimate empirical variograms (actually correlation functions) and then model them using a number of different approaches .. mostly using Matern as basis
  #\\ returns empirical variogram and parameter estimates, and the models themselves
  #\\ expect xy = c(p/lon, p/lat), z= variable
  #\\ varZ is the total variance which needs to be mulitplied to the curve if you want the "true" semivariance
  #\\ parameterization is as in spBayes and gstat
  #\\ covariogram (||x||) = tau^2 * (2^{nu-1} * Gamma(nu) )^{-1} * (phi*||x||)^{nu} * K_{nu}(phi*||x||)
  # -------------------------


  if ( 0 ) {
   # just for debugging / testing ... and example of access method:
   bioLibrary("bio.utilities", "bio.spacetime")
   require(sp)
   data(meuse)
    xy = meuse[, c("x", "y")]
    z = log( meuse$zinc )
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
    hist( out$spBayes$recover$p.theta.samples[,3] ) # phi
    hist( out$spBayes$recover$p.theta.samples[,4] ) # nu

    gr = spacetime.variogram( xy, z, methods="geoR" )
    gs = spacetime.variogram( xy, z, methods="gstat" )
    grf = spacetime.variogram( xy, z, methods="RandomFields" )
    gsp = spacetime.variogram( xy, z, methods="spBayes" )

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

  names(xy) =  c("plon", "plat" )

  drange = sqrt((diff(range(xy$plon)))^2 + (diff(range(xy$plat)))^2)
  drange = min( diff(range(xy$plon)), diff(range(xy$plat)), drange )

  # if max dist not given, make a sensible choice
  if ( is.na(maxdist)) {
    maxdist = drange * 0.5  # default
  } else if ( maxdist=="all") {
    maxdist = drange
  }

  xrange = range( xy$plon, na.rm=TRUE )
  yrange = range( xy$plat, na.rm=TRUE )

  z = z / sd( z, na.rm=TRUE)
  zrange = range( z, na.rm=TRUE )

  difx = diff( xrange)
  dify = diff( yrange)

  dd = max( difx, dify )

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

    model = ~ RMspheric( var=NA, scale=NA) + RMnugget(var=NA) + RMtrend(NA)
    
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

    # kappa  == 1/phi.geoR
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

  if ("LaplacesDemon" %in% methods){
    
    stop("This is too slow to use. It is just to document the approach. dmvn is the culprit .. try sparse matrix methods ...")

    require(LaplacesDemon)
    require(LaplacesDemonCpp)
    
    coords = data.frame( plon=(xy$plon-yrange[1] )/dd,  plat=(xy$plat-xrange[1])/dd )
    
    LD = list(
      N = length(z),
      X = as.matrix(data.frame( intercept=rep(1, length(z)) ) ),
      y = z,  
      D.N = as.matrix(dist( coords, diag=TRUE, upper=TRUE)) #  distance matrix
    )

    LD$parm.names = as.parm.names( list(
      beta=1, 
      sigma=rep(0,2),
      phi=0,
      kappa=1,
      zeta=rep(0,LD$N)
    ) )

    LD$mon.names = c("LP" )
    
    LD$PGF = function(LD) {
      beta = rnorm(1)
      sigma = runif(2,0.1,10)
      phi = runif(1,1,5)
      kappa = 1
      zeta = mvtnorm::rmvnorm(1, rep(0,LD$N), sigma[2]*sigma[2]*exp(-phi*LD$D.N)^kappa, method="chol")
      return(c(zeta, beta, sigma, phi))
    }

    LD$i = list(
      beta = grep("\\<beta\\>", LD$parm.names),
      sigma = grep("\\<sigma\\>", LD$parm.names),
      phi = grep("\\<phi\\>", LD$parm.names),
      kappa = grep("\\<kappa\\>", LD$parm.names),
      zeta = grep("\\<zeta\\>", LD$parm.names)
    )
    
    LD$Model = function(parm, LD){
      
      ### Parameters
      loglik = c()
      loglik[LD$parm.names] = 0

      parm[LD$i$sigma] = sigma = interval(parm[LD$i$sigma], 1, Inf)
      parm[LD$i$phi] = phi = interval(parm[LD$i$phi], 1, 5)
      Sigma = sigma[2]*sigma[2] * exp(-phi * LD$D.N)^parm[LD$i$kappa]
      
      ### Log-Prior
      loglik[LD$i$beta] = sum(dnormv(parm[LD$i$beta], 0, 1000, log=TRUE))
      loglik[LD$i$zeta] = mvtnorm::dmvnorm(parm[LD$i$zeta], rep(0, LD$N), Sigma, log=TRUE)
      loglik[LD$i$sigma] = sum(dhalfcauchy(sigma - 1, 25, log=TRUE))
      loglik[LD$i$phi] = dunif(phi, 1, 5, log=TRUE)
      
      ### Log-Likelihood
      mu = tcrossprod(LD$X, t(parm[LD$i$beta]) ) + parm[LD$i$zeta]

      LL = sum(dnorm(LD$y, mu, sigma[1], log=TRUE))
      
      ### Log-Posterior
      LP = LL + sum(loglik)
      Modelout = list(LP=LP, Dev=-2*LL, Monitor=c(LP), yhat=rnorm(length(mu), mu, sigma[1]), parm=parm)
      return(Modelout)
    }

    Model = compiler::cmpfun(LD$Model) 

    # first mcmc .. faster start-up and convergence to global equil 
    f = LaplacesDemon(Model, Data=LD, Initial.Values=LD$PGF(LD), Iterations=500, Status=1, Thinning=1 )
    
    # now use the currest parm est and try a laplace:
    f = LaplaceApproximation(Model, Data=LD, parm=as.initial.values(f), Iterations=50, Method="LBFGS", CPUs=6 ) # refine it
    
    # f = VariationalBayes(Model, Data=LD, parm=as.initial.values(f), Iterations=500, Samples=20, CPUs=5, Covar=f$Covar ) # refine it again
    # f = LaplacesDemon(Model, Data=LD, Initial.Values=as.initial.values(f), Iterations=1000, Status=100, Thinning=1, Covar=f$Covar )
    # f = VariationalBayes(Model, Data=LD, parm=as.initial.values(f), Iterations=100, Samples=10, CPUs=5, Covar=f$Covar )
    # f = IterativeQuadrature(Model, Data=LD, parm=as.initial.values(f), Iterations=10, Algorithm="AGH",
    #  Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL), Covar=f$Covar )
    
    if (plotdata) {
      Consort(f)
      plot(f)
    }
  }

 
  return(out)
}


