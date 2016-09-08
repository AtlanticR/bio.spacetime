
spacetime.reshape = function(data, locsout, nr, nc, interp.method="kernel.density", trimprobs=c(0.025, 0.975), theta=NULL, nsd=10 ) {
  #\\ reshape statistics to fit the output resolution 
  # trim quaniles in case of extreme values
  
  vq = quantile( data$z, probs=trimprobs, na.rm=TRUE )
  ii = which( data$z < vq[1] )
  if ( length(ii)>0) data$z[ii] = vq[1]
  jj = which( data$z > vq[2] )
  if ( length(jj)>0) data$z[jj] = vq[2]
  # default is base::I = no transform 
  data$z = transform(data$z))
  out = res = NULL
  if (interp.method == "kernel.density") {
    # default :: create a "surface" and reshape to a grid using (gaussian) kernel-based smooth via FFT
    require(fields)
    isurf = fields::interp.surface( surface, loc=locsout  )
    zout = matrix( isurf, nrow=nr, ncol=nc )
    out = fields::image.smooth( zout, theta=theta, xwidth=theta*nsd, ywidth=theta*nsd ) # nsd SD of the normal kernel
  }
  if ( !is.null( out )) res = revTransform( out$z )
  return (res)
}
   
    
     
