
spacetime_interpolate_xy_simple_multiple = function( ip=NULL, p ) {
  #// designed to be called from sapcetime_interpolate
  #// for the sake of speed and parallelization, the kernel density method is written out again 
  #// within the for loop rather than reusing the spacetime_interpolate_simple

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  # spatial interpolation by finest time scale 
  Z0 = matrix( NaN, nrow=p$nplons, ncol=p$nplats)

  for ( iip in ip ) {
    ww = p$runs[ iip, "tiyr_index" ]

    # mean estimates    
    Z = Z0
    P = attach.big.matrix( p$ptr$P )  # predictions
    
    tofill = which( ! is.finite( P[,ww] ) )
    if (length( tofill) == 0 ) next()

    Z[p$Mat2Ploc] = P[,ww]
    Zp = image.smooth( Z, dx=p$pres, dy=p$pres, wght=p$wgts )$z 
    P[,ww][tofill] = Zp[p$Mat2Ploc][ tofill]

    # update and repeat --- check if required ..
    tofill = which( ! is.finite( P[,ww] ) )
    if (length( tofill) > 0 ) {
      Z[p$Mat2Ploc] = P[,ww] # update matrix
      Zp = image.smooth( Z, dx=p$pres, dy=p$pres, wght=p$wgts )$z 
      P[,ww][tofill] = Zp[p$Mat2Ploc][ tofill]
    }

    ## SD estimates
    Z = Z0
    Psd = attach.big.matrix( p$ptr$Psd )   # predictions

    tofill = which( ! is.finite( Psd[,ww] ) )
    if (length( tofill) == 0 ) next()

    Z[p$Mat2Ploc] = Psd[,ww]
    Zp = image.smooth( Z, dx=p$pres, dy=p$pres, wght=p$wgts )$z 
    Psd[,ww][tofill] = Zp[p$Mat2Ploc][ tofill]

    # update and repeat --- check if required ..
    tofill = which( ! is.finite( Psd[,ww] ) )
    if (length( tofill) > 0 ) {
      Z[p$Mat2Ploc] = Psd[,ww] # update matrix
      Zp = image.smooth( Z, dx=p$pres, dy=p$pres, wght=p$wgts )$z 
      Psd[,ww][tofill] = Zp[p$Mat2Ploc][ tofill]
    }
  }

  return( "complete" )

}

