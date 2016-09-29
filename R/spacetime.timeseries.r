
spacetime.timeseries = function( p, YiU, Si, newdata ) {


  if ( p$tsmethod %in% c( "harmonics.1", "harmonics.2", "harmonics.3"  ) ) {
    #\\ simple GAM with spatial weights (inverse distance squared) and harmonics
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
    Sloc = p$ff$Sloc  # statistical output locations
    Yloc = p$ff$Yloc  # statistical output locations

    x = data.frame( p$ff$Y[YiU], p$ff$Ytime[YiU] )
    names(x) = c(p$variables$Y, p$variables$TIME)
    if ( exists("spacetime.link", p) ) x[, p$variables$Y] = p$spacetime.link ( x[, p$variables$Y] ) 

    x$yr = trunc( x[, p$variables$TIME])
    x$plon = Yloc[YiU,1]
    x$plat = Yloc[YiU,2]
    if (exists("COV", p$variables)) {
      x[ p$variables$COV ] = p$ff$Ycov[YiU]
      close(p$ff$Ycov)
    }
    Y_wgt = 1 / (( Sloc[Si,1] - x$plat)**2 + (Sloc[Si,2] - x$plon)**2 )# weight data in space: inverse distance squared
    Y_wgt[ which( Y_wgt < 1e-3 ) ] = 1e-3
    Y_wgt[ which( Y_wgt > 1 ) ] = 1

    close(Yloc)
    close(Ytime)
    close(p$ff$Y)
    close(Sloc)

    x$cos.w  = cos( 2*pi*x$tiyr )
    x$sin.w  = sin( 2*pi*x$tiyr )

    newdata$cos.w  = cos( newdata$tiyr )
    newdata$sin.w  = sin( newdata$tiyr )
    
    newdata$yr = trunc( newdata[,p$variables$TIME] )

    # compute aSiitional harmonics only if required (to try to speed things up a bit)
    if ( p$tsmethod %in% c( "harmonics.2", "harmonics.3"  ) ) {
      x$cos.w2 = cos( 2*x$tiyr )
      x$sin.w2 = sin( 2*x$tiyr )
      newdata$cos.w2 = cos( 2*newdata$tiyr )
      newdata$sin.w2 = sin( 2*newdata$tiyr )
    }
    if ( p$tsmethod %in% c( "harmonics.3"  ) ) {
      x$cos.w3 = cos( 3*x$tiyr )
      x$sin.w3 = sin( 3*x$tiyr )
      newdata$cos.w3 = cos( 3*newdata$tiyr )
      newdata$sin.w3 = sin( 3*newdata$tiyr )
    }
    
    # estimate model parameters
    tsmodel = try( 
      gam( t ~ s(yr) + s(yr, cos.w) + s(yr, sin.w) + s(cos.w) + s(sin.w), 
        data=x, weights=Y_wgt, optimizer=c("outer","bfgs")  ) ) 

    if ( class(tsmodel) %in% "try-error" ) next()
   
    out = try( predict( tsmodel, newdata=newdata, type="response", se.fit=T ) ) 

    if ( ! "try-error" %in% class( out ) ) {
      newdata$fit = as.vector(out$fit)
      newdata$sd = as.vector(out$se.fit) # this is correct: se.fit== stdev of the mean fit: eg:  https://stat.ethz.ch/pipermail/r-help/2005-July/075856.html
      if (exists("spacetime.invlink", p)) {
        newdata$fit =  p$spacetime.invlink( newdata$fit )
        newdata$sd  =  p$spacetime.invlink( newdata$sd )
      }
      return(newdata)  
    } else {
      return( NULL )
    }

  }


  if (p$tsmethod=="") {
    newdata = NULL

    return(newdata) 
  }
  

}
