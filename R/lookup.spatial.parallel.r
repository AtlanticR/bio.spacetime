
  lookup.spatial.parallel = function( X, Y, init.files, distance.threshold.km=1, ix=1, iy=2, iz=3, cltype="SOCK", clusters="localhost" ) {
    nid = nrow( X )
    cl = makeCluster( spec=clusters, type=cltype)
    ssplt = lapply( clusterSplit(cl, 1:nid), function(i) i )   # subset data into lists
    pp.out = clusterApplyLB( cl, ssplt, lookup.spatial, X=X, Y=Y, init.files=init.files, distance.threshold.km=1, ix=ix, iy=iy, iz=iz )
    stopCluster(cl)
    X = NULL
    for (m in 1:length(pp.out))  X = rbind( X, pp.out[[m]] )
    return ( X )
  }


