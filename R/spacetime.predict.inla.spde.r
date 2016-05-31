
spacetime.predict.inla.spde = function( MESH, RES, p, stack_index=NULL, locs_new=NULL, task=NULL ) {
  
  out = NULL

  if ( "predictions.direct" %in% task ) {
    # precomputed ... slow and expensive in RAM/CPU, just extract from tag indices
    xmean = RES$summary.fitted.values[ stack_index, "mean"]
    xsd = RES$summary.fitted.values[ stack_index, "sd"]
    if (exists("spacetime.invlink", p)) {
      xmean =  p$spacetime.invlink( xmean )
      xsd =  p$spacetime.invlink( xsd )
    }
    out = data.frame(xmean=xmean, xsd=xsd) 
    return(out)
  }


  if ( "predictions.projected" %in% task ) {
    #\\ note this method only works with simple additive models 
    #\\ when smoothes are involved, it becomes very complicated and direct estimation is probably faster/easier
    pG = inla.mesh.projector( MESH, loc=as.matrix( locs_new) )
    posterior.samples = inla.posterior.sample(n=p$inla.nsamples, RES)
    rnm = rownames(posterior.samples[[1]]$latent )  
    posterior = sapply( posterior.samples, p$spacetime.posterior.extract, rnm=rnm )
    if (exists("spacetime.invlink", p)) posterior = p$spacetime.invlink( posterior )   # return to user scale
    rm(posterior.samples); gc()
    # robustify the predictions by trimming extreme values .. will have minimal effect upon mean
    # but variance estimates should be useful/more stable as the tails are sometimes quite long 
    for (ii in 1:nrow(posterior )) {
      qnt = quantile( posterior[ii,], probs=p$predict.quantiles, na.rm=TRUE ) 
      toolow = which( posterior[ii,] < qnt[1] )
      toohigh = which (posterior[ii,] > qnt[2] )
      if (length( toolow) > 0 ) posterior[ii,toolow] = qnt[1]
      if (length( toohigh) > 0 ) posterior[ii,toohigh] = qnt[2]
    }
    xmean = c( inla.mesh.project( pG, field=apply( posterior, 1, mean, na.rm=TRUE )  ))
    xsd   = c( inla.mesh.project( pG, field=apply( posterior, 1, sd, na.rm=TRUE )  ))
    out = data.frame(xmean=xmean, xsd=xsd) 
    return( out )
  }
}


