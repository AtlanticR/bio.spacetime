
spacetime__bayesx = function( p, x, pa ) {
  #\\ this is the core engine of spacetime .. localised space-time modelling interpolation and prediction .. using bayesx 
   
  # EG: 
  #  logzinc ~  sx( x,y, bs="kr")  # "kr" is perhaps overly smooth
  #  logzinc ~  sx( x,y, bs="te")  # more detail .. "te" is preferred

  if ( !exists( "bayesx.method", p) ) p$bayesx.method="MCMC"  # slightly more smoothing than the REML method

  hmod = try( bayesx( p$spacetime_engine_modelformula, data=x, method=p$bayesx.method ) )

  if ( "try-error" %in% class(hmod) ) return( NULL )

  px = predict(hmod)
  ss = summary(lm( px ~ x[, p$variables$Y ], na.action="na.omit" ))
  if (ss$r.squared < p$spacetime_rsquared_threshold ) return(NULL)
    
  out = try( predict( hmod, newdata=pa, type="response" ) ) 


  # plot( hmod, term = "sx(x,y)", image=TRUE)
  # summary(hmod)
  # lattice::levelplot( out ~ plon+plat, data=k, aspect="iso" )

  if ( "try-error" %in% class( out ) ) return( NULL )

  pa$mean = out
  pa$sd = 1 # no option right now to estim posterior prediction errors .. may be possible with type="terms" but would be slow to simulate  and do not know how to do it yet .. please fix this ..

  spacetime_stats = list( sdTotal=sd(x[,p$variable$Y], na.rm=T), rsquared=ss$r.squared, ndata=nrow(x) ) # must be same order as p$statsvars .. pseudo rsquared for logistic .. for poisson {1- logLik(mod) / logLik(mod_saturated)} might be better
  
  # lattice::levelplot( mean ~ plon + plat, data=pa[pa$tiyr==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  
  return( list( predictions=pa, spacetime_stats=spacetime_stats ) )  
}
