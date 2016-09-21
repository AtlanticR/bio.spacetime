
spacetime.covariance.spatial = function( ip=NULL, p ) {
  #\\ estimate spatially localized spatial covariance parameters

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns
  
  #---------------------
  # data for modelling
  # dependent vars
  Y = p$ff$Y  # for read only
  hasdata = 1:nrow(Y)
  bad = which( !is.finite( Y[]))
  if (length(bad)> 0 ) hasdata[bad] = NA

  # data locations
  Yloc = p$ff$Yloc # read only
  bad = which( !is.finite( rowSums(Yloc[])))
  if (length(bad)> 0 ) hasdata[bad] = NA

  hasdata = na.omit(hasdata)
  Yloc_good = Yloc[hasdata,]

  #-----------------
  # row, col indices for statistical outputs
  Sloc = p$ff$Sloc    # statistical output locations, read only
  rcS = data.frame( cbind( 
      Srow = (Sloc[,1]-p$plons[1])/p$pres + 1,  
      Scol = (Sloc[,2]-p$plats[1])/p$pres + 1))

  # main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    dd = p$runs[ iip, "jj" ]
    # print (dd)
    S = p$ff$S  # inside loop to update the results immediately ; read-write
    if ( is.nan( S[dd,1] ) ) next()
    if ( !is.na( S[dd,1] ) ) next()
    S[dd,1] = NaN  # flag: if a run fails, do not revisit, over-written below if successful
    # choose a distance <= p$dist.max where n is within range of reasonable limits to permit a numerical solution
    ppp = NULL
    ppp = try( point.in.block( Sloc[dd,], Yloc_good, 
      dist.max=p$dist.max, dist.min=p$dist.min, n.min=p$n.min, n.max=p$n.max,
      upsampling=p$upsampling, downsampling=p$downsampling, resize=TRUE ) )
    if( is.null(ppp)) next()
    if (class( ppp ) %in% "try-error" ) next()
    dist.cur = ppp$dist.to.nmax
    j = hasdata[ppp$indices]
    rm(ppp)
    ndata = length(j) # number of data locations
    if (ndata < p$n.min) next()
    xy = as.data.frame(  )
    z = Y[j]

    print( "Computing variogram" )
    res = NULL
    res = spacetime.variogram( Yloc[j,], z, methods=p$variogram.engine )

    if (!is.null(res)) {
      if (exists(p$variogram.engine, res) ) {
        print( "Saving summary statisitics" )
        # save statistics last as this is an indicator of completion of all tasks .. restarts would be broken otherwise
        S[dd,1] = res$varZ
        S[dd,2] = res[[p$variogram.engine]]$varSpatial
        S[dd,3] = res[[p$variogram.engine]]$varObs
        S[dd,4] = res[[p$variogram.engine]]$range
        S[dd,5] = res[[p$variogram.engine]]$phi
        S[dd,6] = res[[p$variogram.engine]]$nu
    }}


    if(0) {
      x11();
      levelplot( ( S[,4] ) ~ plon + plat, data=Sloc[,], aspect="iso",
        labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }
    
    close(S)

  }  # end for loop

  close(Sloc)
  close(Yloc)
  close(Y)

  return( "complete" )

}


