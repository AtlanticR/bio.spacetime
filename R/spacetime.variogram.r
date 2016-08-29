
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
    maxdist = NA
    edge=c(1/3, 1)
    nbreaks = 15

    nc_max = 5  # max number of iterations

    out = list()
    out$varZ = var( z, na.rm=TRUE )  # this is the scaling factor for semivariance .. diving by sd, below reduces numerical floating point issues
    out$meanZ = mean(z, na.rm=TRUE)
    out$minX = min( xy[,1], na.rm=TRUE )
    out$minY = min( xy[,2], na.rm=TRUE )

    #scaling xyz helps stabilize and speed up solutions
    z = (z - out$meanZ )/ sqrt( out$varZ ) # (centered and scaled by sd to avoid floating point issues)
    zrange = range( z, na.rm=TRUE )   

    names(xy) =  c("plon", "plat" ) # arbitrary
    xr = range( xy$plon, na.rm=TRUE )
    yr = range( xy$plat, na.rm=TRUE )
    drange = min( diff( xr), diff( yr)  )

    # if max dist not given, make a sensible choice
    if ( is.na(maxdist)) {
      maxdist = drange * 0.1  # default
    } else if ( maxdist=="all") {
      maxdist = drange
    }
    out$drange = drange
    out$maxdist = maxdist

    # positions and distances are scaled to max dist ..
    xy$plon = ( xy$plon - out$minX ) / maxdist
    xy$plat = ( xy$plat - out$minY ) / maxdist
  }

  if (0) {
    # tests
    gr = spacetime.variogram( xy, z, methods="geoR" )
    gs = spacetime.variogram( xy, z, methods="gstat" )
    grf = spacetime.variogram( xy, z, methods="RandomFields" )
    gsp = spacetime.variogram( xy, z, methods="spBayes" )
    ginla = spacetime.variogram( xy, z, methods="inla" )

    # tests:
    out = gsp
    nd = nrow(out$spBayes$recover$p.theta.samples)
    rr = rep(NA, nd )
    for (i in 1:nd) rr[i] = geoR::practicalRange("matern", phi=1/out$spBayes$recover$p.theta.samples[i,3], kappa=out$spBayes$recover$p.theta.samples[i,4] )
    hist(rr)  # range estimate

    hist( out$spBayes$recover$p.theta.samples[,1] ) #"sigma.sq"
    hist( out$spBayes$recover$p.theta.samples[,2] ) # "tau.sq"
    hist( out$spBayes$recover$p.theta.samples[,3] ) # 1/phi
    hist( out$spBayes$recover$p.theta.samples[,4] ) # nu

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

  # ------------------------

  if ("gstat" %in% methods){
    #\\ covariogram (||x||) = tau^2 * (2^{nu-1} * Gamma(nu) )^{-1} * (phi*||x||)^{nu} * K_{nu}(phi*||x||)
    #\\ gstat::kappa == spBayes::nu
    #\\ gstat::range == spBayes::phi {aka, "scale parameter"}

    require(gstat)
    require(sp)

    vrange = 0.5 # starting est of range
    distx = vrange * 0.9 ## back it up a bit to enter smoothly into the loop
    nc = 0
    while ( distx < vrange ) {
      nc = nc  + 1
      distx = distx * 1.25 # gradually increase distx until solution found
      vEm = try( variogram( z~1, locations=~plon+plat, data=xy, cutoff=distx, width=distx/nbreaks, cressie=TRUE ) ) # empirical variogram
      if  ("try-error" %in% vEm) return(NULL)
      vMod0 = vgm(psill=0.75, model="Mat", range=distx, nugget=0.25, kappa=1 ) # starting model parameters
      #vMod0 = vgm("Mat")
      vFitgs =  try( fit.variogram( vEm, vMod0, fit.kappa =TRUE, fit.sills=TRUE, fit.ranges=TRUE ) ) ## gstat's kappa is the Bessel function's "nu" smoothness parameter
      vrange = min(1, geoR::practicalRange("matern", phi=vFitgs$range[2], kappa=vFitgs$kappa[2]  ) )
      if (nc > nc_max ) break()
    }

    if  ("try-error" %in% vFitgs) return(NULL)
    vEm$dist = vEm$dist * out$maxdist
    vEm$gamma = vEm$gamma * out$varZ
    vFitgs$psill = vFitgs$psill * out$varZ
    vFitgs$range[2] = out$maxdist* vFitgs$range[2]

    out$gstat = list( fit=vFitgs, vgm=vEm, range=NA, nu=vFitgs$kappa[2], phi=vFitgs$range[2],
        varSpatial=vFitgs$psill[2], varObs=vFitgs$psill[1]  )  # gstat::"range" == range parameter == phi
    
    out$gstat$range = geoR::practicalRange("matern", phi=out$gstat$phi, kappa=out$gstat$nu  )

    if (plotdata) {
      x11()
      plot(vEm, model=vFitgs, add=T)
      x11()
      plot( gamma ~ dist, data=out$gstat$vgm, xlim=c(0,maxdist), 
           ylim=c(0,max(out$gstat$vgm$gamma)*1.1), col="blue", pch=20 )
      abline( h=out$gstat$varSpatial + out$gstat$varObs ) 
      abline( h=out$gstat$varObs )
      abline( v=out$gstat$range )
      x = seq( 0, maxdist, length.out=100 )
      acor = geoR::matern( x, phi=out$gstat$phi, kappa=out$gstat$nu  )
      acov = out$gstat$varObs + out$gstat$varSpatial * (1- acor)
      lines( acov~x , col="red" )
    }
    return(out)
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
    vrange = 0.5
    distx = vrange * 0.9 ## back it up a bit to enter smoothly into the loop
    nc = 0
    while ( distx  < vrange ) {
      nc = nc + 1
      distx = distx * 1.25
      vEm = try( variog( coords=xy, data=z, uvec=nbreaks, max.dist=distx ) )
      if  ("try-error" %in% vEm) return(NULL)
      vMod = try( variofit( vEm, nugget=0.5, kappa=1, cov.model="matern", ini.cov.pars=c(0.5, distx/4) ,
        fix.kappa=FALSE, fix.nugget=FALSE, max.dist=distx, weights="cressie" ) )
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
  
    out$geoR = list( fit=vMod, vgm=vEm, model=vMod, range=vMod$practicalRange*out$maxdist,
              varSpatial= vMod$cov.pars[1]*out$varZ, varObs=vMod$nugget*out$varZ, 
              nu=vMod$kappa,  phi=vMod$cov.pars[2]*out$maxdist )

    if (plotdata) {
      # not rescaled ...
      x11()
      plot(vEm)
      lines(vMod)
      x11()
      plot( out$geoR$vgm )
      x11()
      plot( out$geoR$vgm$v*out$varZ ~ c(out$geoR$vgm$u*out$maxdist), pch=20 , 
           xlim=c(0,maxdist), ylim=c(0, out$varZ*1.25) )
      abline( h=out$geoR$varSpatial + out$geoR$varObs)  
      abline( h=out$geoR$varObs )
      abline( v=out$geoR$range )
      x = seq( 0, max(out$geoR$vgm$u*out$maxdist), length.out=100 )
      acor = geoR::matern( x, phi=out$geoR$phi, kappa=out$geoR$nu  )
      acov = out$geoR$varObs +  out$geoR$varSpatial * (1-acor)  ## geoR is 1/2 of gstat and RandomFields gamma's
      lines( acov ~ x , col="orange" )
    }

    return(out)

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
   # RandomFields:  Cov(h) = v * Orig(A*h/s) ; s=scale, h=dist, A=aniso, v=variance, Original model (data scale)

   # where nu > 0 and K_nu is the modified Bessel function of second kind and distance r >= 0 between two pointsd
   # The Matern covariance model is given by: C(h) = v * phi(A*h/s).
   #  Cov(r) = 2^{1- nu} Gamma(nu)^{-1} (sqrt{2nu} r)^nu K_nu(sqrt{2nu} r)
   # "phi" = sqrt{2nu}/?? NOt clear ...
   
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
              varSpatial=oo$param["value", "matern.var"]*out$varZ,
              varObs=oo$param["value", "nugget.var"]*out$varZ,
              phi=(out$maxdist* oo$param["value", "matern.s"] )/(sqrt(oo$param["value", "matern.nu"]*2) ), 
              nu=oo$param["value", "matern.nu"], # RF::nu == geoR:: kappa (bessel smoothness param)
              error=NA )

    out$RandomFields$range = geoR::practicalRange("matern", phi=out$RandomFields$phi, kappa=out$RandomFields$nu  )

    if (plotdata) {
      x11()
      py = as.vector(out$RandomFields$vgm@emp.vario) *out$varZ 
      px = out$RandomFields$vgm@centers*out$maxdist
      plot(  py ~ px, pch=20, ylim=c(0,out$varZ*1.25) )
      abline( h=out$RandomFields$varSpatial + out$RandomFields$varObs  )
      abline( h=out$RandomFields$varObs )
      abline( v=out$RandomFields$range )

      x = seq( 0, max(px ), length.out=100 )
      acor = geoR::matern( x, phi=out$RandomFields$phi, kappa=out$RandomFields$nu  )
      acov = out$RandomFields$varObs  +  out$RandomFields$varSpatial*(1- acor)
      lines( acov~x , col="red" )

      # compare with:
      x11()
      plot(o)
    }
    return(out)

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

    burn.in <- 0.2*n.samples

    ##recover beta and spatial random effects
    m.1 <- spRecover(model, start=burn.in )

    u = apply(m.1$p.theta.recover.samples, 2, mean)
    u["phi"] = out$maxdist/u["phi"]
    u["sigma.sq"] = u["sigma.sq"]*out$varZ
    u["tau.sq"] = u["tau.sq"]*out$varZ

    vrange = geoR::practicalRange("matern", phi=u["phi"], kappa=u["nu"]  )

    out$spBayes = list( model=model, recover=m.1,
      range=vrange, varSpatial=u["sigma.sq"], varObs=u["tau.sq"], 
      phi=u["phi"], nu=u["nu"] )  # output using geoR nomenclature

    if (plotdata) {
      x11()
      x = seq( 0, vrange* 2, length.out=100 )
      acor = geoR::matern( x, phi=u["phi"], kappa=u["nu"] )
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

    return(out)

  }


  # -------------------------


  if ("inla" %in% methods){
    require(INLA)
    require(lattice)

    inla.setOption(scale.model.default = TRUE)  # better numerical performance of IGMRF models and less dependnence upon hyperpriors

    locs0  = as.matrix( xy )
    xy$b0 = 1  # intercept for inla

    vRange = maxdist * 0.1

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
        # control.compute=list(dic=TRUE),
        control.results=list(return.marginals.random=TRUE ),
        # control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        # control.fixed = list(expand.factor.strategy='inla') ,
        control.predictor=list(A=inla.stack.A(Z), compute=TRUE, link=1 ) ,
        # control.inla = list( h=1e-4, tolerance=1e-10),
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

    out$inla = list(summary=inla.summary, 
      mesh=MESH, res=RES, range.inla90=inla.summary[["range","mean"]]*out$maxdist,
      varSpatial=inla.summary[["spatial error","mean"]]*out$varZ, 
      varObs=inla.summary[["observation error","mean"]]*out$varZ,
      phi = out$maxdist/inla.summary[["kappa","mean"]] , nu=alpha-1, error=NA )

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
  
    return(out)

  }


  # -------------------------

  if ("BayesX" %in% methods){
    library("R2BayesX")
    # fixes nu=1.5
    # phi = max(distance) / const, such that Corr(distance=const) = 0.001; 
    # i.e. range at distance where covar ~0.999 .. but not sure how to recover the correct phi/range from this ...

    fm1 <- bayesx( z ~ sx(plon, plat, bs="kr" ), family="gaussian", method="REML", data =xy )
    out$BayesX = list( fit=fm1, range=NA, jagsmodel=fm1, 
        varSpatial=fm1$smooth.hyp[,"Variance"]*out$varZ, 
        varObs=fm1$variance*out$varZ, 
        nu=1.5, phi=NA )
    # out$BayesX$range = geoR::practicalRange("matern", phi=out$BayesX$phi, kappa=out$BayesX$nu  )
    return(fm1)
  }


  # -------------------------

  if ("jags" %in% methods){
    require(rjags)
    require(jagsUI)
    # assume nu = 1 (due to identifiability issues)

    print( "Slow ... 7.5 min for meuse test data")

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
       n.iter=1000,
       n.chains=3,
       n.burnin=100,
       n.thin=5,
       parallel=TRUE,
       DIC=FALSE)
   
   if (0) {
     summary(fit)
     plot(fit)
     gelman.plot(fit$samples)
    # geweke.plot(fit$samples)
    #update(fit, n.iter=2000, n.thin=20 )
      acf( fit$sims.list$phi)
      acf( fit$sims.list$sigmasq)
      acf( fit$sims.list$tausq)

# JAGS output for model '/tmp/RtmpXoFpO7/file539863b2591e', generated by jagsUI.
# Estimates based on 3 chains of 2000 iterations,
# burn-in = 200 iterations and thin rate = 10,
# yielding 540 total samples from the joint posterior. 
# MCMC ran in parallel for 7.461 minutes at time 2016-08-10 18:05:54.

#          mean    sd  2.5%   50% 97.5% overlap0 f  Rhat n.eff
# phi     1.252 0.597 0.513 1.094 2.784    FALSE 1 1.014   238
# sigmasq 0.441 0.084 0.292 0.436 0.621    FALSE 1 1.008   356
# tausq   0.145 0.096 0.030 0.121 0.374    FALSE 1 1.013   279
# $phi
# [1] 222.3888

# $sigmasq
# [1] 0.08303909

# $tausq
# [1] 0.02724344

# $range
# [1] 889.2267

    }
    print (summary(fit))

    out$jags = list(
      fit = fit, model=jagsmodel,
      phi = out$maxdist / fit$summary["phi", "mean"],
      sigmasq = fit$summary["sigmasq", "mean"]*out$varZ,
      tausq = fit$summary["tausq", "mean"]*out$varZ
    )
    out$jags$range = geoR::practicalRange("matern", phi=out$jags$phi, kappa=1  )
    return(out)
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
      eps = 1e-6,
      N = length(z),  # required for LaplacesDemon
      DIST=as.matrix(dist( xy, diag=TRUE, upper=TRUE)), # distance matrix between knots
      y=z  
    )
    Data$mon.names = c( "LP", paste0("yhat[",1:Data$N,"]" ) )
    Data$parm.names = as.parm.names(list(tausq=0, sigmasq=0, phi=0, nu=0 ))
    Data$pos = list(
      tausq = grep("tausq", Data$parm.names),
      sigmasq = grep("sigmasq", Data$parm.names),
      phi = grep("phi", Data$parm.names),
      nu = grep("nu", Data$parm.names)
    )
    Data$PGF = function(Data) {
      #initial values .. get them near the center of mass
      tausq = rgamma (1, 1, 5) # 0 to 1.5 range
      sigmasq = rgamma (1, 1, 5)
      phi = rgamma (1, 1, 1)  # 0 to 500 range
      nu = runif(1, 0.5, 4)
      return( c( tausq, sigmasq, phi, nu ))
    }
    Data$PGF  = compiler::cmpfun(Data$PGF)
    Data$Model = function(parm, Data) {
      tausq = parm[Data$pos$tausq] = LaplacesDemonCpp::interval_random(parm[Data$pos$tausq], Data$eps, 1, 0.01 )
      sigmasq = parm[Data$pos$sigmasq]= LaplacesDemonCpp::interval_random(parm[Data$pos$sigmasq], Data$eps, 1, 0.01 )
      phi = parm[Data$pos$phi]= LaplacesDemonCpp::interval_random(parm[Data$pos$phi], Data$eps, Inf, 0.01 )
      nu = parm[Data$pos$nu] = LaplacesDemonCpp::interval_random(parm[Data$pos$nu], 0.1, 4.0, 0.01 )
      # corSpatial = exp(-Data$DIST/phi)^nu   ## spatial correlation .. exponential
      corSpatial = geoR::matern( Data$DIST, phi=1/phi, kappa=nu )   ## spatial correlation .. matern
      # corSpatial = zapsmall(corSpatial)
      #uphi <- Data$DIST/phi
      #corSpatial = (((2^(-(nu-1)))/gamma(nu)) * (uphi^nu) * besselK(x=phi, nu=nu))
      #diag(corSpatial) = 1
      #if (any(!is.finite(corSpatial))) browser()
      # browser()
      eSp = rmvn( 1, rep(0, Data$N), sigmasq*corSpatial )# psill
      eObs = rnorm( Data$N, 0, sqrt(tausq) ) # nugget error
      tausq.prior = dgamma(tausq, 1, 1, log=TRUE) # 0-1.55 range
      sigmasq.prior = dgamma(sigmasq, 1, 1, log=TRUE)
      phi.prior = dgamma(phi, 1, 1, log=TRUE)
      nu.prior = dnorm(nu, 1, 0.01, log=TRUE)
      yhat = eObs + eSp # local iid error + spatial error
      LL = sum(dnorm(Data$y, yhat, sqrt(sigmasq+tausq), log=TRUE)) ## Log Likelihood
      LP = sum(LL, sigmasq.prior, tausq.prior, phi.prior, nu.prior) ### Log-Posterior
      Modelout = list(LP=LP, Dev=-2*LL, Monitor=c(LP, yhat), yhat=yhat, parm=parm)
      return(Modelout)
    }
    Data$Model = compiler::cmpfun(Data$Model) #  byte-compiling for more speed .. use RCPP if you want more speed
    
    parm0=Data$PGF(Data)
  
    f = LaplaceApproximation(Data$Model, Data=Data, parm=parm0, Method="HAR", Iterations=10000, CovEst="Hessian", sir=TRUE, Interval=1e-8, Samples=5000, Stop.Tolerance=1e-8 )

    f = LaplaceApproximation(Data$Model, Data=Data, parm=parm0, Method="BFGS", Iterations=10000, CovEst="Hessian", sir=TRUE, Interval=1e-6, Samples=5000, Stop.Tolerance=1e-6 )

    f = LaplaceApproximation(Data$Model, Data=Data, parm=parm0, Method="CG", Iterations=5000, CovEst="Hessian", sir=TRUE, Interval=1e-9 ) 

    mu = f$Summary1[,1]
    f0 = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), 
      Iterations=10000, Thinning=10, Status=1000, Algorithm="IM", Specs=list(mu=mu), 
      Covar=f$Covar, CPUs=8 )


    mu = apply(f0$Posterior1, 2, mean)
    f0 = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), 
      Iterations=10000, Thinning=10, Status=1000, Algorithm="IM", Specs=list(mu=mu), 
      Covar=f$Covar, CPUs=8 )

    if (plotdata) {
      
    parm0 = as.initial.values(f)

      f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=parm0, CPUs=8)
      f = LaplaceApproximation(Data$Model, Data=Data, parm=parm0, Method="SPG", Iterations=500, CovEst="Hessian", sir=FALSE )    

      Consort(f)
      plot(f, Data=Data)
       m = f$Summary2[grep( "\\<yhat\\>", rownames( f$Summary2 ) ),]
      # m = f$Summary2[grep( "muSpatial", rownames( f$Summary2 ) ),]
      plot( Data$y ~ m[, "Mean"]  )

      Consort(f0)
      plot(f0, Data=Data)
    }

    out$LaplacesDemon = list( fit=f, vgm=NA, model=Data$Model, range=NA,
      varSpatial=f$Summary2["sigmasq", "Mean"] *out$varZ, 
      varObs=f$Summary2["tausq", "Mean"]*out$varZ, 
      nu=f$Summary2["nu", "Mean"],  
      phi=out$maxdist/f$Summary2["phi", "Mean"] 
    )
    out$LaplacesDemon$range = geoR::practicalRange("matern", phi=out$LaplacesDemon$phi, kappa=out$LaplacesDemon$nu)
 
    out$LaplacesDemon

    if (plotdata) {
      x11()
      x = seq( 0,  out$LaplacesDemon$range * 1.25, length.out=100 )
      svar =  out$LaplacesDemon$varObs + out$LaplacesDemon$varSpatial * (1-geoR::matern( x, phi=out$LaplacesDemon$phi, kappa=out$LaplacesDemon$nu  ))
      plot( svar~x, type="l", ylim=c(0, max(svar)) )
      abline( h=out$LaplacesDemon$varObs + out$LaplacesDemon$varSpatial )
      abline( h=out$LaplacesDemon$varObs )
      abline( v=out$LaplacesDemon$range, col="red"  )
    }
  
  return(out)
  
  }

}


