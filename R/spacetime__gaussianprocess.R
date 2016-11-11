
spacetime__gaussianprocess = function( p, x, pa, timeslices=1, method="fields.Krig" ) {
  #\\ this is the core engine of spacetime .. localised space (no-time) modelling interpolation 
  # \ as a gaussian process
  #\\ note: time is not being modelled and treated independently 
  #\\      .. you had better have enough data in each time slice

  if (!exists("fields.covariancefunction", p)) p$fields.covariancefunction="Matern"  # note: "Rad.cov" is TPS
  if (!exists("nu", p)) p$nu=2  # note: this is the smoothness or shape parameter (fix at 2 if not calculated or given )   
  if (!exists("tps_p", p)) p$tps_p=2 # power in the radial basis function (with a log term added for 
     # even dimensions). 
     # If m is the degree of derivative in penalty then p=2m-d 
     # where d is the dimension of x. p must be greater than 0. 
   
  x_r = range(x[,p$variables$LOCS[1]])
  x_c = range(x[,p$variables$LOCS[2]])

  x_nr = diff(x_r)/p$pres + 1
  x_nc = diff(x_c)/p$pres + 1

  x_plons = seq( x_r[1], x_r[2], length.out=x_nr )
  x_plats = seq( x_c[1], x_c[2], length.out=x_nc )

  x_locs = expand.grid( x_plons, x_plats ) # final output grid
  attr( x_locs , "out.attrs") = NULL
  names( x_locs ) = p$variables$LOCS

  x$pred = NA

  # locations of the new (output) coord system .. smaller than the data range of x
  pa_r = range(pa[,p$variables$LOCS[1]])
  pa_c = range(pa[,p$variables$LOCS[2]])
  
  pa_nr = diff(pa_r)/p$pres + 1
  pa_nc = diff(pa_c)/p$pres + 1

  pa_plons = seq( pa_r[1], pa_r[2], length.out=pa_nr )
  pa_plats = seq( pa_c[1], pa_c[2], length.out=pa_nc )
  
  pa_locs = expand.grid( pa_plons, pa_plats ) # final output grid
  attr( pa_locs , "out.attrs") = NULL
  names( pa_locs ) = p$variables$LOCS

  x$mean = NA

  pa$mean = NA
  pa$sd = NA

  for ( ti in timeslices ) {
    
    if ( exists("TIME", p$variables) ) {
      xi = which( x[ , p$variables$TIME ] == ti )
    } else {
      xi = 1:nrow(x) # all data as p$nt==1
    }

    xy = x[xi, p$variables$LOCS]
    z = x[xi, p$variables$Y]
    fsp = try( MLESpatialProcess.fast(xy, z, cov.function = "stationary.cov",  cov.args = list(Covariance=p$fields.covariancefunction, smoothness=p$nu ) ) )

    if (inherits(vFitgs, "try-error") )  next()
    if ( fsp$converge != 0 ) next()

    lambda.MLE<- fsp$pars[3]/fsp$pars[1]  # ratio nugget / sill variance
    fspmodel <- Krig( xy, z, Covariance=p$fields.covariancefunction, theta=fsp$pars[2], smoothness=p$nu, lambda= lambda.MLE, p=p$tps_p )

    x$mean[xi] = predict(fspmodel, x=x[xi, p$variables$LOCS] )
    ss = lm( x$mean[xi] ~ x[xi,p$variables$Y], na.action=na.omit)
    if ( "try-error" %in% class( ss ) ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$spacetime_rsquared_threshold ) next()

    pa_i = which( pa[, p$variables$TIME]==ti)
    pa$mean[pa_i] = predict(fspmodel, x=pa[pa_, p$variables$LOCS] )
    pa$sd[pa_i]   = predict(fspmodel, x=pa[pa_, p$variables$LOCS] )

    if ( 0 ){
      # debugging plots
      surface(fsp2)
      fsp.p<- predictSurface(fsp2, lambda= lambda.MLE, nx=200, ny=200, )
      surface(fsp.p, type="I")
      fsp.p2<- predictSurfaceSE(fsp2,)
      surface(fsp.p2, type="C")
    }
 
  }

  # plot(pred ~ z , x)
  ss = lm( x$mean ~ x[,p$variables$Y], na.action=na.omit)
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$spacetime_rsquared_threshold ) return(NULL)

  spacetime_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, spacetime_stats=spacetime_stats ) )  
}

