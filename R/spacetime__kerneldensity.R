
spacetime__kerneldensity = function( p, x, pa ) {
  #\\ this is the core engine of spacetime .. localised space (no-time) modelling interpolation 
  #\\ note: time is not being modelled and treated independently 
  #\\      .. you had better have enough data in each time slice

  x_r = range(x[,p$variables$LOCS[1]])
  x_c = range(x[,p$variables$LOCS[2]])

  x_nr = diff(x_r)/p$pres + 1
  x_nc = diff(x_c)/p$pres + 1

  x_plons = seq( x_r[1], x_r[2], length.out=x_nr )
  x_plats = seq( x_c[1], x_c[2], length.out=x_nc )

  x_locs = expand.grid( x_plons, x_plats ) # final output grid
  attr( x_locs , "out.attrs") = NULL
  names( x_locs ) = p$variables$LOCS

  x$mean = NA

  pa$mean = NA
  pa$sd = NA

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

  for ( ti in p$timeslices ) {
    
    if ( exists("TIME", p$variables) ) {
      xi = which( x[ , p$variables$TIME ] == ti )
    } else {
      xi = 1:nrow(x) # all data as p$nt==1
    }

    # map of row, col indices of input data in the new (output) coordinate system
    l2M = cbind( ( x[xi,p$variables$LOCS[1]]-x_r[1])/p$pres + 1, 
                  (x[xi,p$variables$LOCS[2]]-x_c[1])/p$pres + 1 )
   
    # matrix representation of the output surface
    M = matrix( NA, nrow=x_nr, ncol=x_nc) 
    M[l2M] = x[xi,p$variables$Y] # fill with data in correct locations

    stats = rep( NA, nrow( pa_locs) )  # output data
       
    Z = try( fields::image.smooth( M, dx=p$pres, dy=p$pres, theta=p$theta) )
    if ( "try-error" %in% class(Z) ) next()

    # match prediction to input data 
    x_id = cbind( ( x[xi,p$variables$LOCS[1]]-x_r[1])/p$pres + 1, 
                   (x[xi,p$variables$LOCS[2]]-x_c[1])/p$pres + 1 )
    x$mean[xi] = Z$z[x_id]

    ss = lm( x$mean[xi] ~ x[xi,p$variables$Y], na.action=na.omit)
    if ( "try-error" %in% class( ss ) ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$spacetime_rsquared_threshold ) next()

    pa_i = which( pa[, p$variables$TIME]==ti)
    Z_i = cbind( ( pa[pa_i,p$variables$LOCS[1]]-x_r[1])/p$pres + 1, 
                  (pa[pa_i,p$variables$LOCS[2]]-x_c[1])/p$pres + 1 )

    # make sure predictions exist .. kernel density can stop prediction beyond a given range if the xwidth/ywidth options are not used and/or the kernel distance (theta) is small 
    if ( any( Z_i<1) ) next()  
    if ( any( Z_i[,1] > x_nr) ) next()
    if ( any( Z_i[,2] > x_nc) ) next()

    pa$mean[pa_i] = Z$z[Z_i]
    pa$sd[pa_i] = 1
  }

  # plot(mean ~ z , x)
  ss = lm( x$mean ~ x[,p$variables$Y], na.action=na.omit )
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$spacetime_rsquared_threshold ) return(NULL)

  spacetime_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=rsquared, ndata=nrow(x) ) # must be same order as p$statsvars
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, spacetime_stats=spacetime_stats ) )  
}

