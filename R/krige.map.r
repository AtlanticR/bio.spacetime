
krige.map = function( p=NULL) {
  if (!p$do.parallel) {
    krige.map.core( p=p)
  } else  {
    cl = makeCluster( spec=p$clusters, type=p$cltype)
    ssplt = lapply( clusterSplit(cl, 1:p$nruns), function(i) i )   # subset data into lists
    clusterApplyLB( cl, ssplt, krige.map.core, p=p)
    stopCluster(cl)
  }
}

