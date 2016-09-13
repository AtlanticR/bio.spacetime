
spacetime.covariance.spatial = function( ip=NULL, p ) {
  #\\ estimate spatially localized spatial covariance parameters

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns
  p = spacetime.db( p=p, DS="filenames" )

  #---------------------
  # data for modelling
  # dependent vars
  Y = h5file( p$ptr$Y ) 
  hasdata = 1:length(Y)
  bad = which( !is.finite( Y[]))
  if (length(bad)> 0 ) hasdata[bad] = NA

  # data locations
  Yloc = h5file( p$ptr$Yloc )
  bad = which( !is.finite( rowSums(Yloc[])))
  if (length(bad)> 0 ) hasdata[bad] = NA

  hasdata = na.omit(hasdata)
  Yloc_good = Yloc[hasdata,]

  #-----------------
  # row, col indices for statistical outputs
  Sloc = h5file( p$ptr$Sloc )  # statistical output locations
  rcS = data.frame( cbind( 
      Srow = (Sloc[,1]-p$plons[1])/p$pres + 1,  
      Scol = (Sloc[,2]-p$plats[1])/p$pres + 1))

  # main loop over each output location in S (stats output locations)
  S = h5file( p$ptr$S )
  for ( iip in ip ) {
    dd = p$runs[ iip, "jj" ]
    focal = t(Sloc[dd,])
    # print (dd)

    if ( is.nan( S[dd,1] ) ) next()
    if ( !is.na( S[dd,1] ) ) next()

    S[dd,1] = NaN   # this is a flag such that if a run fails (e.g. in mesh generation), it does not get revisited
    # .. it gets over-written below if successful
    # choose a distance <= p$dist.max where n is within range of reasonable limits to permit a numerical solution
    # slow ... need to find a faster solution

    ppp = NULL
    ppp = try( point.in.block( focal[1,c(1,2)], Yloc_good, dist.max=p$dist.max, dist.min=p$dist.min, n.min=p$n.min, n.max=p$n.max,
      upsampling=p$upsampling, downsampling=p$downsampling, resize=TRUE ) )

    if( is.null(ppp)) next()
    if (class( ppp ) %in% "try-error" ) next()
    dist.cur = ppp$dist.to.nmax

    j = hasdata[ppp$indices]
    rm(ppp )
    ndata = length(j) # number of data locations
    # print ( paste( "Ndata :", ndata ) )
    if (ndata < p$n.min) next()

    xy = as.data.frame( Yloc[j,] )
    z = Y[j]

    print( "Computing variogram" )
    res = NULL
    res = spacetime.variogram( xy, z, methods=p$variogram.engine )

    if (!is.null(res)) {
      if (exists(p$variogram.engine, res) ) {
        print( "Saving summary statisitics" )
        # save statistics last as this is an indicator of completion of all tasks .. restarts would be broken otherwise
        S[dd,1] = res$varZ
        S[dd,2] = res[[p$variogram.engine]]$varSpatial
        S[dd,3] = res[[p$variogram.engine]]$varObs
        S[dd,4] = res[[p$variogram.engine]]$range
        S[dd,5] = res[[p$variogram.engine]]$phi
        S[dd,6] = res[[p$variogram.engine]]$kappa
    }}

    if(0) {
      pps = expand.grid( plon=p$plons, plat=p$plats)
      zz = which(pps$plon > min(pa$plon) & pps$plon < max(pa$plon) & pps$plat < max(pa$plat) & pps$plat > min(pa$plat) )
      x11();
      levelplot( ( S[zz,"range"] ) ~ plon + plat, pps[zz,], aspect="iso",
        labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
    }


  }  # end for loop
  h5close(S)
  h5close(Sloc)
  h5close(Yloc)
  h5close(Y)

  return( "complete" )

}


