
krige.map = function( p=NULL, init.files=NULL ) {
  if (!p$do.parallel) {
    krige.map.core( p=p, init.files=init.files )
  } else  {
    cl = makeCluster( spec=p$clusters, type=p$cltype)
    ssplt = lapply( clusterSplit(cl, 1:p$nruns), function(i) i )   # subset data into lists
    clusterApplyLB( cl, ssplt, krige.map.core, p=p, init.files=init.files )
    stopCluster(cl)
  }
}

