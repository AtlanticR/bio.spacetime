 
  # ---------------------------------------------------------------------
  # Wrapping function to access GMT (Generic Mapping Tools)
  #     
  #     Warning: GMT calls must be understood by your system: 
  #     i.e., GMT program must be in the file path 
  #     "Gawk" is used in some instances ..  this can eventually be implemented directly in R
  
  gmt.map.variables = function( U, params, variables, plottimes, basedir, conversions="ps2png", delta=1, init.files=NULL, db="snowcrab", cltype="SOCK") {
    require(parallel)
    nid = length( variables )
    if (!params$do.parallel) {
      gmt.map.variables.core(U=U, params=params, variables=variables, plottimes=plottimes, basedir=basedir, conversions=conversions, delta=delta, init.files=init.files, db=db)
    } else if (params$do.parallel) {
      cl = makeCluster( spec=params$clusters, type=cltype)
      ssplt = lapply( clusterSplit(cl, 1:nid), function(i) i )   # subset data into lists
      clusterApplyLB( cl, ssplt, gmt.map.variables.core, 
        U=U, params=params, variables=variables, plottimes=plottimes, basedir=basedir, conversions=conversions, delta=delta, init.files=init.files, db=db )
      stopCluster(cl)
    }
    return ("Completed mapping")
  }


