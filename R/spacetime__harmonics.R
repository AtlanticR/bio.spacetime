
spacetime__harmonics = function( p, YiU, Si, pa ) {
   #\\ this is the core engine of spacetime .. localised space-time modelling
   #\\ the simplest form is gam/kernel-based 

    #\\ simple GAM with spatial weights (inverse distance squared) and ts harmonics
    # forcing sinusoid as a seasonal component: 
    # to add an offset to a trig function (b) must add cos to a sin function
    # y ~ a + c*sin(x+b)
    # y ~ a + c*sin(b)*cos(x) + c*cos(b)*sin(x)  
    #   .. as C*sin(x+b) = C*( cos(b) * sin(x) + sin(b) * cos(x) )
    # y ~ b0 + b1*x1 + b2*x2
    # where: 
    #   a = b0
    #   c^2 = b1^2 + b2^2 = c^2*(sin^2(b) + cos^2(b))
    #   c = sqrt(b1^2 + b2^2)
    #   b1/b2 = tan(b)  
    #   b = arctan(b1/b2)
    
    Sloc =  p$ptr$Sloc 
    Yloc =  p$ptr$Yloc 
    Y =  p$ptr$Y 

    x = data.frame( Y[YiU] )
    names(x) = p$variables$Y
    if ( exists("spacetime.link", p) ) x[, p$variables$Y] = p$spacetime.link ( x[, p$variables$Y] ) 
    x$plon = Yloc[YiU,1]
    x$plat = Yloc[YiU,2]
    x$Y_wgt = 1 / (( Sloc[Si,1] - x$plat)**2 + (Sloc[Si,2] - x$plon)**2 )# weight data in space: inverse distance squared
    x$Y_wgt[ which( x$Y_wgt < 1e-3 ) ] = 1e-3
    x$Y_wgt[ which( x$Y_wgt > 1 ) ] = 1
    
    if (exists("COV", p$variables)) {
      Ycov = ( p$ptr$Ycov )
      for (i in 1:length(p$variables$COV )) x[, p$variables$COV[i] ] = Ycov[YiU,i]
    }
    
    if (exists("TIME", p$variables)) {
      Ytime = ( p$ptr$Ytime )
      x[, p$variables$TIME ] = Ytime[YiU,] 
      x$yr = trunc( x[, p$variables$TIME])
      x$cos.w  = cos( 2*pi*x$tiyr )
      x$sin.w  = sin( 2*pi*x$tiyr )
      pa$yr = trunc( pa[,p$variables$TIME] )
      pa$cos.w  = cos( pa$tiyr )
      pa$sin.w  = sin( pa$tiyr )
      # compute adiitional harmonics only if required (to try to speed things up a bit)
      if ( p$spacetime_engine %in% c( "harmonics.2", "harmonics.3"  ) ) {
        x$cos.w2 = cos( 2*x$tiyr )
        x$sin.w2 = sin( 2*x$tiyr )
        pa$cos.w2 = cos( 2*pa$tiyr )
        pa$sin.w2 = sin( 2*pa$tiyr )
      }
      if ( p$spacetime_engine %in% c( "harmonics.3"  ) ) {
        x$cos.w3 = cos( 3*x$tiyr )
        x$sin.w3 = sin( 3*x$tiyr )
        pa$cos.w3 = cos( 3*pa$tiyr )
        pa$sin.w3 = sin( 3*pa$tiyr )
      }
    }

    # estimate model parameters
    hmod = try( 
      gam( p$spacetime_engine_modelformula, data=x, weights=Y_wgt, optimizer=c("outer","bfgs")  ) ) 

    if ( "try-error" %in% class(hmod) ) next()
    
    out = try( predict( hmod, newdata=pa, type="response", se.fit=T ) ) 

    if ( "try-error" %in% class( out ) ) return( NULL )

    pa$mean = as.vector(out$fit)
    pa$sd = as.vector(out$se.fit) # this is correct: se.fit== stdev of the mean fit: eg:  https://stat.ethz.ch/pipermail/r-help/2005-July/075856.html

    if (exists("spacetime.invlink", p)) {
      pa$mean = p$spacetime.invlink( pa$mean )
      pa$sd  =  p$spacetime.invlink( pa$sd )
    }

    if (exists( "quantile_bounds", p)) {
      tq = quantile( x[,p$variables$Y], probs=p$quantile_bounds, na.rm=TRUE  )
      bad = which( pa$mean < tq[1] | pa$mean > tq[2]  )
      if (length( bad) > 0) {
        pa$mean[ bad] = NA
        pa$sd[ bad] = NA
      }
    }

    ss = summary(hmod)
    spacetime_stats = list( sdTotal=sd(Y[], na.rm=T), rsquared=ss$r.sq, ndata=ss$n ) # must be same order as p$statsvars

    return( list( predictions=pa, spacetime_stats=spacetime_stats ) )  

}
