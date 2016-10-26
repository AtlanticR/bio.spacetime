
spacetime__harmonics = function( p, x, pa ) {
  #\\ this is the core engine of spacetime .. localised space-time modelling interpolation and prediction
  #\\ simple GAM with spatial weights (inverse distance squared) and ts harmonics 
  
  # estimate model parameters
  hmod = try( 
    gam( p$spacetime_engine_modelformula, data=x, weights=Y_wgt, optimizer=c("outer","bfgs")  ) ) 

  if ( "try-error" %in% class(hmod) ) next()
  
  out = try( predict( hmod, newdata=pa, type="response", se.fit=T ) ) 

  if ( "try-error" %in% class( out ) ) return( NULL )

  pa$mean = as.vector(out$fit)
  pa$sd = as.vector(out$se.fit) # this is correct: se.fit== stdev of the mean fit: eg:  https://stat.ethz.ch/pipermail/r-help/2005-July/075856.html

  ss = summary(hmod)
  spacetime_stats = list( sdTotal=sd(Y[], na.rm=T), rsquared=ss$r.sq, ndata=ss$n ) # must be same order as p$statsvars

  return( list( predictions=pa, spacetime_stats=spacetime_stats ) )  

}
