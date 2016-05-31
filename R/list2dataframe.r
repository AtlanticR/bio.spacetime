
  list2dataframe = function( x, split=",") {
    y = strsplit(as.character( x ), split=split )
    nr = length( y )
    z = unlist( y )
    nc = length (z) / nr
    o = as.data.frame( matrix( z, nrow=nr, ncol=nc, byrow=T) )
    return( o )
  }



