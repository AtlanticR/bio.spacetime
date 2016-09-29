
point_in_block = function( Sloc, Si, Yloc, Yi, dist.max, dist.min, n.min=30, n.max=1000, upsampling=c(1.1, 1.2), downsampling=c(0.9, 0.8), resize=FALSE ) {

  #\\ data within a specified range/distance/number of points using ff data backend

  U = NULL
  dlon = abs( Sloc[Si,1] - Yloc[Yi,1] ) 
  dlat = abs( Sloc[Si,2] - Yloc[Yi,2] ) 
  U = which( dlon  <= dist.max  & dlat <= dist.max ) # faster to take a block 
  ndat = length(U)
  dist.cur = dist.max
  
  if ( !resize )  returnreturn(list(U=U, dist=dist.cur))

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
