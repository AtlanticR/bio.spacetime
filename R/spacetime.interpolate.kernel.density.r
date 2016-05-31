

spacetime.interpolate.kernel.density = function( x, y, z, locsout, nxout, nyout, theta, xwidth, ywidth ) {
  # create a "surface" and interpolate to a different grid using (gaussian) kernel-based smooth 
  require(fields)
  isurf = fields::interp.surface( list( x=x, y=y, z=z), loc=locsout  )
  zout = matrix( isurf, nrow=nxout, ncol=nyout )
  res = fields::image.smooth( zout, theta=theta, xwidth=xwidth, ywidth=ywidth ) 
  return( res ) 
}

  

