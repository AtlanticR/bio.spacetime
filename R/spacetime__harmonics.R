
spacetime__harmonics = function( p, YiU, Si, newdata ) {
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
    
    Sloc = attach.big.matrix( p$ptr$Sloc )
    Yloc = attach.big.matrix( p$ptr$Yloc )
    Y = attach.big.matrix( p$ptr$Y )

    x = data.frame( Y[YiU] )
    names(x) = p$variables$Y
    if ( exists("spacetime.link", p) ) x[, p$variables$Y] = p$spacetime.link ( x[, p$variables$Y] ) 
    x$plon = Yloc[YiU,1]
    x$plat = Yloc[YiU,2]
    Y_wgt = 1 / (( Sloc[Si,1] - x$plat)**2 + (Sloc[Si,2] - x$plon)**2 )# weight data in space: inverse distance squared
    Y_wgt[ which( Y_wgt < 1e-3 ) ] = 1e-3
    Y_wgt[ which( Y_wgt > 1 ) ] = 1
    
    if (exists("COV", p$variables)) {
      Ycov = attach.big.matrix( p$ptr$Ycov )
      for (i in 1:length(p$variables$COV )) x[, p$variables$COV[i] ] = Ycov[YiU,i]
    }
    
    if (exists("TIME", p$variables)) {
      Ytime = attach.big.matrix( p$ptr$Ytime )
      x[, p$variables$TIME ] = Ytime[YiU,] 
      names(x) = c(p$variables$Y, p$variables$TIME)
      x$yr = trunc( x[, p$variables$TIME])
      x$cos.w  = cos( 2*pi*x$tiyr )
      x$sin.w  = sin( 2*pi*x$tiyr )
      newdata$yr = trunc( newdata[,p$variables$TIME] )
      newdata$cos.w  = cos( newdata$tiyr )
      newdata$sin.w  = sin( newdata$tiyr )
      # compute adiitional harmonics only if required (to try to speed things up a bit)
      if ( p$spacetime_engine %in% c( "harmonics.2", "harmonics.3"  ) ) {
        x$cos.w2 = cos( 2*x$tiyr )
        x$sin.w2 = sin( 2*x$tiyr )
        newdata$cos.w2 = cos( 2*newdata$tiyr )
        newdata$sin.w2 = sin( 2*newdata$tiyr )
      }
      if ( p$spacetime_engine %in% c( "harmonics.3"  ) ) {
        x$cos.w3 = cos( 3*x$tiyr )
        x$sin.w3 = sin( 3*x$tiyr )
        newdata$cos.w3 = cos( 3*newdata$tiyr )
        newdata$sin.w3 = sin( 3*newdata$tiyr )
      }
    }

 
    # estimate model parameters
    tsmodel = try( 
      gam( p$modelformula, data=x, weights=Y_wgt, optimizer=c("outer","bfgs")  ) ) 

    if ( "try-error" %in% class(tsmodel) ) next()
    
    out = try( predict( tsmodel, newdata=newdata, type="response", se.fit=T ) ) 

    if ( "try-error" %in% class( out ) ) return( NULL )

    newdata$mean = as.vector(out$fit)
    newdata$sd = as.vector(out$se.fit) # this is correct: se.fit== stdev of the mean fit: eg:  https://stat.ethz.ch/pipermail/r-help/2005-July/075856.html

    if (exists("spacetime.invlink", p)) {
      newdata$mean = p$spacetime.invlink( newdata$mean )
      newdata$sd  =  p$spacetime.invlink( newdata$sd )
    }

    if (exists( "Y_bounds", p)) {
      bad = which( newdata$mean < p$Y_bounds[1] | newdata$mean > p$Y_bounds[2]  )
      if (length( bad) > 0) {
        newdata$mean[ bad] = NA
        newdata$sd[ bad] = NA
      }
    }

    if (exists( "quantile_bounds", p)) {
      tq = quantile( Y, probs=p$quantile_bounds, na.rm=TRUE  )
      bad = which( newdata$mean < tq[1] | newdata$mean > tq[2]  )
      if (length( bad) > 0) {
        newdata$mean[ bad] = NA
        newdata$sd[ bad] = NA
      }
    }

    ss = summary(tsmodel)
    spacetime_stats = list( sdTotal=sd(Y[], na.rm=T), rsquared=ss$r.sq, ndata=ss$n ) # must be same order as p$statsvars

    return( list( predictions=newdata, spacetime_stats=spacetime_stats ) )  

}
