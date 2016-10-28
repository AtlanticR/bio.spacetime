
# NOTE:: not finishe porting fully ... only designed for xy data right no .. no time .. needs more testing

  spacetime__inla = function( p, x, pa ) {
    #\\ generic spatial and space-time interpolator using inla
    #\\ parameter and data requirements can be seen in bathymetry\src\bathymetry.r
  
    # the following parameters are for inside and outside ... do not make them exact multiples as this seems to make things hang ..
    if ( !exists("inla.mesh.max.edge", p))  p$inla.mesh.max.edge = c(  0.025,   0.04 )    # proportion of 2*p$dist.max or equivalent: c(inside,outside) -- must be positive valued
    if ( !exists("inla.mesh.offset", p))  p$inla.mesh.offset   = c( - 0.025,  - 0.05 )   # how much to extend inside and outside of boundary: proportion of dist.max .. neg val = proportion
    if ( !exists("inla.mesh.cutoff", p)) p$inla.mesh.cutoff   = c( - 0.05,   - 0.5 )    ## min distance allowed between points: proportion of dist.max ; neg val = proportion

    if ( !exists("inla.mesh.hull.radius", p)) p$inla.mesh.hull.radius = c( -0.04, - 0.08 ) ## radius of boundary finding algorythm ; neg val = proportion

    if ( !exists("inla.mesh.hull.resolution", p)) p$inla.mesh.hull.resolution = 125  ## resolution for discretization to find boundary

    if ( !exists("spacetime.noise", p)) p$spacetime.noise = 0.001  # add a little noise to coordinates to prevent a race condition

    if ( !exists("inla.alpha", p)) p$inla.alpha = 2 # bessel function curviness
    if ( !exists("inla.nsamples", p)) p$inla.nsamples = 5000 # posterior similations 
    if ( !exists("predict.in.one.go", p)) p$predict.in.one.go = FALSE # use false, one go is very very slow and a resource expensive method
    if ( !exists("predict.quantiles", p)) p$predict.quantiles = c(0.025, 0.975 )  # posterior predictions robustified by trimming extreme values 
    
    if ( !exists("debug.file", p)) p$debug.file = file.path( bio.workdirectory, "inla.debug.out" )

    # priors
    # kappa0 = sqrt(8)/p$expected.range
    # tau0 = 1/(sqrt(4*pi)*kappa0* p$expected.sigma)

    #-----------------
    # row, col indices


      locs_noise = runif( nrow(x)*2, min=-p$pres*p$spacetime.noise, max=p$pres*p$spacetime.noise ) # add  noise  to prevent a race condition

      # also sending direct distances rather than proportion seems to cause issues..
      MESH = spacetime_mesh( locs=x[,p$variables$LOC] + locs_noise,
        lengthscale=p$spacetime_distance_prediction*2,
        max.edge=p$inla.mesh.max.edge * p$spacetime_distance_prediction*2,
        bnd.offset=p$inla.mesh.offset,
        cutoff=p$inla.mesh.cutoff,
        convex=p$inla.mesh.hull.radius,
        resolution=p$inla.mesh.hull.resolution )

      rm(locs_noise)

      if (0)  plot( MESH )  # or ... spacetime_plot( p=p, "mesh", MESH=MESH )

      SPDE = inla.spde2.matern( MESH,
        alpha=p$inla.alpha #, # alpha is the Bessel smoothness factor .. 1(?) gives exponential correlation function
        # B.tau=matrix(c(log(tau0),-1,+1),nrow=1,ncol=3),
        # B.kappa=matrix(c(log(kappa0),0,-1),nrow=1,ncol=3),
        # theta.prior.mean=c(0,0), # thetas are the internal representations of prior offsets for tau and kappa (i.e.,range and sigma)
        # theta.prior.prec=c(0.1, 0.1) # precision of thetas
      )

      # effects .. a list of two elements: first is for SPDE and second is for covariates
      obs_index = inla.spde.make.index(name="spatial.field", SPDE$n.spde)
      obs_eff = list()
      obs_eff[["spde"]] = c( obs_index, list(intercept=1) )
      if ( exists( "COV", p$variables) ) {
        covar = as.data.frame( x[, p$variables$COV ] ) 
        colnames( covar ) = p$variables$COV
        obs_eff[["covar"]] = as.list(covar)
        obs_A = list( inla.spde.make.A( mesh=MESH, loc=x[, p$variables$LOC ] ), 1 )
      } else {
        obs_A = list( inla.spde.make.A( mesh=MESH, loc=x[, p$variables$LOC ] ) ) # no effects
      }

      obs_ydata = list()
      obs_ydata[[ p$variables$Y ]] = x[, p$variables$Y ]
      
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

      # prediction stack:: check for covariates

        preds_locs = as.matrix( pa[, p$variables$LOC ] )
        preds_index = inla.spde.make.index( name="spatial.field", SPDE$n.spde)
        preds_eff = list()
        preds_eff[["spde"]] = c( preds_index, list(intercept=1) )
        if ( exists( "COV", p$variables) ) {
          pcovars = as.data.frame(pa[,p$variables$COV])
          colnames( pcovars ) = p$variables$COV
          preds_eff[["covar"]] = as.list( pcovars )
          preds_A = list( inla.spde.make.A(MESH, loc=preds_locs ), 1)
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

      
      if (!exists("spacetime_engine_modelformula", p) )  p$spacetime_engine_modelformula = formula( z ~ -1 + intercept + f( spatial.field, model=SPDE ) ) # SPDE is the spatial covariance model .. defined in 

      RES = NULL
      RES = spacetime_inla_call( FM=p$spacetime_engine_modelformula, DATA=DATA, SPDE=SPDE, FAMILY=as.character(p$spacetime_family) )
  
      if (is.null(RES))  return(NULL)

      # inla.spde2.matern creates files to disk that are not cleaned up:
      spdetmpfn = SPDE$f$spde2.prefix
      fns = list.files( dirname( spdetmpfn ), all.files=TRUE, full.names=TRUE, recursive=TRUE, include.dirs=TRUE )
      oo = grep( basename(spdetmpfn), fns )
      if(length(oo)>0) file.remove( sort(fns[oo], decreasing=TRUE) )
      rm ( DATA, oo, fns, spdetmpfn); gc()

      # -----------
      # predictions
      
        preds = NULL
        if (!(exists, "spacetime_inla_prediction", p)) p$spacetime_inla_prediction="direct"

        if ( p$spacetime_inla_prediction=="direct" ) {
          # precomputed ... slow and expensive in RAM/CPU, just extract from tag indices
          pa$mean = RES$summary.fitted.values[ preds_stack_index, "mean"]
          pa$sd = RES$summary.fitted.values[ preds_stack_index, "sd"]
        }

        if ( p$spacetime_inla_prediction=="projected") {
          #\\ note this method only works with simple additive models 
          #\\ when smoothes are involved, it becomes very complicated and direct estimation is probably faster/easier
          pG = inla.mesh.projector( MESH, loc=as.matrix( pa[,p$variables$LOC]) )
          posterior.samples = inla.posterior.sample(n=p$inla.nsamples, RES)
          rnm = rownames(posterior.samples[[1]]$latent )  
          posterior = sapply( posterior.samples, p$spacetime.posterior.extract, rnm=rnm )
          rm(posterior.samples); 
          # robustify the predictions by trimming extreme values .. will have minimal effect upon mean
          # but variance estimates should be useful/more stable as the tails are sometimes quite long 
          for (ii in 1:nrow(posterior )) {
            qnt = quantile( posterior[ii,], probs=p$predict.quantiles, na.rm=TRUE ) 
            toolow = which( posterior[ii,] < qnt[1] )
            toohigh = which (posterior[ii,] > qnt[2] )
            if (length( toolow) > 0 ) posterior[ii,toolow] = qnt[1]
            if (length( toohigh) > 0 ) posterior[ii,toohigh] = qnt[2]
          }
          pa$mean = c( inla.mesh.project( pG, field=apply( posterior, 1, mean, na.rm=TRUE )  ))
          pa$sd   = c( inla.mesh.project( pG, field=apply( posterior, 1, sd, na.rm=TRUE )  ))
        }

        if (0) {
          x11(); levelplot( mean ~ plon+plat, pa, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
          x11(); levelplot( sd   ~ plon+plat, pa, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        }

      rm(MESH)

  # bathymetry.figures( DS="predictions", p=p ) # to debug
    inla.summary = spacetime_summary_inla_spde2 ( RES, SPDE )
    # save statistics last as this is an indicator of completion of all tasks .. restarts would be broken otherwise
    stats = list()
    stats$sdTotal=sd(x[,p$variable$Y], na.rm=T)
    stats$rsquared=NA
    stats$ndata=nrow(x)
    stats$varSpatial = inla.summary["spatial error", "mode"]
    stats$varObs = inla.summary["observation error", "mode"]
    stats$range = inla.summary["range", "mode"]
    stats$range.sd  = inla.summary["range", "sd"]

    return (list(predictions=pa, spacetime_stats=stats))
  }


