
  lookup.spatial = function( id=NULL, X, Y, init.files, distance.threshold.km=1, ix=1, iy=2, iz=3 ) {
    
    # this can be made faster by using plon and plat and a direct search for distance
    # ie. rectangular blocks but at the cost of warped space at the extremes ... should explore
    
    for (i in init.files) source( i )
    
    coords = c( iy, ix ) ## iy=lon and iy=lat
    ### iz  is column index for variable of interest

    if (is.null (id) )  id = c(1:nrow(X))
    id = as.numeric( id) 

    X$lookup.mean = NA
    X$lookup.var = NA
    X$lookup.nw = NA
    X$lookup.n = NA

    # compute variogram parameters this way too ...     
    for ( i in id ) {
      print(i)
      
      # approximation appropriate for scotian shelf 
      # initial screen based upon 1 nm ~ 1.85 km, or 1 km ~ 0.54 nm
      # 60 nm/deg lon at equator
      # 1 degree lon ~ 80 km or 0.0125 deg lon / km,  .. use a 
      lon.multiplier = 0.0125 * 1.1  # the additional multiplier 1.1 adds a small buffer
      # 1 degree lat ~ 110 km or 0.009 deg lat / km   .. 0.010
      lat.multiplier = 0.009 * 1.1

      dsx = c(-1,1) * distance.threshold.km *  lon.multiplier + X[i,ix] 
      dsy = c(-1,1) * distance.threshold.km * lat.multiplier + X[i,iy] 

      Y.sample = Y [ which( 
        Y[,ix] >= dsx[1] & Y[,ix] <= dsx[2] &
        Y[,iy] >= dsy[1] & Y[,iy] <= dsy[2] 
      ), ]

      d = geosphere::distCosine( X[i, coords], Y.sample[, coords] ) / 1000
      near = NULL
      near = which(d <= distance.threshold.km) # within 1 km of the point
      nn =  length(near) 
      if ( nn>= 1 ) {
        vals = Y.sample[near,iz]
        weights = 1/sqrt(d[near])
        X$lookup.mean[i] = wtd.mean(x=vals, w=weights, na.rm=T)
        X$lookup.var[i] = wtd.var(x=vals, w=weights, na.rm=T)
        X$lookup.nw[i] = sum( weights, na.rm=T)  
        X$lookup.n[i] = nn  
      } 
    }
    
    out = X[id,] # these will be added together if they are in segments
     
    return(out)
  }


