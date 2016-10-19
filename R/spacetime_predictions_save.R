
spacetime_predictions_save = function( p, pred ) {
  # update SD estimates of predictions with those from other locations via the
  # incremental  method ("online algorithm") of mean estimation after Knuth ;
  # see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
  # update means: inverse-variance weighting   
  # see https://en.wikipedia.org/wiki/Inverse-variance_weighting
 
  npred = nrow(pred)

  P = switch( p$storage.backend, 
    bigmemory.ram=attach.big.matrix(p$ptr$P), 
    bigmemory.filebacked=attach.big.matrix(p$ptr$P), 
    ff=p$ptr$P )
  Pn = switch( p$storage.backend, 
    bigmemory.ram=attach.big.matrix(p$ptr$Pn), 
    bigmemory.filebacked=attach.big.matrix(p$ptr$Pn), 
    ff=p$ptr$Pn )
  Psd = switch( p$storage.backend, 
    bigmemory.ram=attach.big.matrix(p$ptr$Psd), 
    bigmemory.filebacked=attach.big.matrix(p$ptr$Psd), 
    ff=p$ptr$Psd )

  if ( ! exists("TIME", p$variables) ) {

    u = which( is.finite( P[pred$i] ) )  # these have data already .. update
    if ( length( u ) > 0 ) {
      ui = pred$i[u]  # locations of P to modify
      Pn[ui] = Pn[ui] + 1 # update counts
      stdev_update =  Psd[ui] + ( pred$sd[u] -  Psd[ui] ) / Pn[ui]
      means_update = ( P[ui] / Psd[ui]^2 + pred$mean[u] / pred$sd[u]^2 ) / ( Psd[ui]^(-2) + pred$sd[u]^(-2) )
      mm = which(is.finite( means_update + stdev_update ))
      if( length(mm)> 0) {
        iumm = ui[mm]
        Psd[iumm] = stdev_update[mm]
        P  [iumm] = means_update[mm]
      }
    }

    # first time # no data yet
    v = setdiff(1:npred, u)         
    if ( length(v) > 0 ) {
      vi = pred$i[v]
      Pn [vi] = 1
      P  [vi] = pred$mean[v]
      Psd[vi] = pred$sd[v]
    }

  }


  if ( exists("TIME", p$variables) ) {

    u = which( is.finite( P[pred$i,1] ) )  # these have data already .. update
    nu = length( u ) 
    if ( nu > 0 ) {
      ui = pred$i[u]  # locations of P to modify
      nc = ncol(P)
      add.ff(Pn, 1, ui, 1:nc ) # same as Pn[ui,] = Pn[ui]+1 but 2X faster
      stdev_update =  Psd[ui,] + ( pred$sd[u] -  Psd[ui,] ) / Pn[ui,]
      means_update = ( P[ui,] / Psd[ui,]^2 + pred$mean[u] / pred$sd[u]^2 ) / 
        ( Psd[ui,]^(-2) + pred$sd[u]^(-2) )
      mm = which( is.finite( rowSums(means_update + stdev_update )))  # created when preds go outside quantile bounds .. this removes all data from a given location rather than the space-time .. severe but likely due to a poor prediction and so remove all (it is also faster this way as few manipulations)
      if( length(mm)> 0) {
        iumm = ui[mm] 
        Psd[iumm] = stdev_update[mm]
        P  [iumm] = means_update[mm]
      }
    }

    # do this as a second pass in case NA's were introduced by the update .. unlikely , but just in case
    v = setdiff(1:npred, u) 
    nv = length(v)          # no data yet
    if ( nv > 0 ) {
      vi = pred$i[v]
      Pn [vi,] = 1
      P  [vi,] = pred$mean[v]
      Psd[vi,] = pred$sd[v]
    }
  
  }
  
  close(P)
  close(Pn )
  close(Psd )

  invisible()

}
