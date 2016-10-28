
point_in_block = function( Sloc, Si, Yloc, Yi, Y, dist.max, dist.min, n.min=30, n.max=1000, upsampling=c(1.1, 1.2), downsampling=c(0.9, 0.8), resize=FALSE, method="simple" ) {

  #\\ data within a specified range/distance/number of points using ff data backend

  U = NULL
  dlon = abs( Sloc[Si,1] - Yloc[Yi,1] ) 
  dlat = abs( Sloc[Si,2] - Yloc[Yi,2] ) 
  U = which( dlon  <= dist.max  & dlat <= dist.max ) # faster to take a block 
  ndat = length(U)
  dist.cur = dist.max
  
  if ( !resize )  return(list(U=U, dist=dist.cur))

  if (method=="variogram" & ndat > n.min ) {
    dat = cbind( Y[Yi], Yloc[Yi,] )
    names(dat) = c("y", "plon", "plat" )
    varY = var(Y[Yi], na.rm=TRUE)
    vEm = try( variogram( y~1, locations=~plon+plat, data=dat, cutoff=dist.max ) ) # empirical variogram
    if  ("try-error" %in% vEm) return(NULL)
    vMod0 = vgm(psill=0.5*varY, model="Sph", range=dist.max*0.75, nugget=0.5*varY ) # starting model parameters
      #vMod0 = vgm("Mat")
    vFitgs =  try( fit.variogram( vEm, vMod0, fit.sills=TRUE, fit.ranges=TRUE ) ) 
    vrange = max( min(1, geoR::practicalRange("sph", phi=vFitgs$range[2] ) ), dist.max)
    U = which( dlon < vrange & dlat < vrange )# faster to take a block 
    return(list(U=U, dist=vrange))
  }

  
  if (method =="simple" | ndat < n.min ) {
    # find data nearest S[Si,] and with sufficient data
    if ( ndat < n.min )  {
      for ( usamp in upsampling )  {
        dist.cur = dist.max * usamp
        U = which( dlon < dist.cur & dlat < dist.cur ) # faster to take a block 
        ndat = length(U)
        if ( ndat >= n.min ) {
          if (ndat >= n.max) U = U[ .Internal( sample( length(U), n.max, replace=FALSE, prob=NULL)) ] 
          return(list(U=U, dist=dist.cur))
        }
      }
    } else if ( ndat >= n.min ) {
      if ( ndat <= n.max * 1.5 ) { # if close to n.max, subsample quickly and return(list(U=U, dist=dist.cur))
        if ( ndat > n.max) { 
          U = U[ .Internal( sample(  length(U), n.max, replace=FALSE, prob=NULL)) ] 
          return(list(U=U, dist=dist.cur))
        } else {
          for ( dsamp in downsampling )  {
            dist.cur = dist.max * dsamp
            U = which( dlon < dist.cur & dlat < dist.cur )# faster to take a block 
            ndat = length(U)
            if ( ndat <= n.max ) return(list(U=U, dist=dist.cur))
            if ( dist.cur <= dist.min ) {
              # reached lower limit in distance, taking a subsample instead
              U = which( dlon < dist.min & dlat < dist.min ) # faster to take a block 
              U = U[ .Internal( sample( length(U), n.max, replace=FALSE, prob=NULL)) ]
              return(list(U=U, dist=dist.cur))
            }
          }
        }
      }
    }
    return(list(U=U, dist=dist.cur))
  }

}
