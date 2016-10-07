
spacetime_interpolate_xy = function( ip=NULL, p ) {
  #\\ estimate spatially localized spatial covariance parameters

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns
  
  #---------------------
  # data for modelling
  # dependent vars
  Y = p$ff$Y  # for read only
  Yi = 1:nrow(Y)
  bad = which( !is.finite( Y[]))
  if (length(bad)> 0 ) Yi[bad] = NA

  # data locations
  Yloc = p$ff$Yloc # read only
  bad = which( !is.finite( rowSums(Yloc[])))
  if (length(bad)> 0 ) Yi[bad] = NA
  Yi = na.omit(Yi)

  #-----------------
  # row, col indices for statistical outputs
  Sloc = p$ff$Sloc    # statistical output locations, read only
  rcS = data.frame( cbind( 
      Srow = (Sloc[,1]-p$plons[1])/p$pres + 1,  
      Scol = (Sloc[,2]-p$plats[1])/p$pres + 1))

  # main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    Si = p$runs[ iip, "locs" ]
    # print (Si)
    S = p$ff$S  # inside loop to update the results immediately ; read-write
    if ( is.infinite( S[Si,1] ) ) next()
    if ( !is.nan( S[Si,1] ) ) next()
    S[Si,1] = Inf  # flag: if a run fails, do not revisit, over-written below if successful
    # choose a distance <= p$dist.max where n is within range of reasonable limits to permit a numerical solution
   
    # find data withing a given distance / number 
    pib = point_in_block( Sloc=Sloc, Si=Si, Yloc=Yloc, Yi=Yi, 
      dist.max=p$dist.max, dist.min=p$dist.min, 
      n.min=p$n.min, n.max=p$n.max, 
      upsampling=p$upsampling, downsampling=p$downsampling, resize=TRUE ) 
    if ( is.null(pib)) {
      next()
    } else {
      dist.cur = pib$dist
      U = pib$U
      rm(pib); gc()
    }
    ndata = length(U)
    if ((ndata < p$n.min) | (ndata > p$n.max) ) next()
    YiU = Yi[U]

    print( "Computing variogram" )
    res = NULL
    res = spacetime_variogram( Yloc[YiU,], Y[YiU], methods=p$variogram.engine )

   # --- predictions
   # do localised predictions here 
    if (p$spacetime_engine=="kernel.density") pred = spacetime_interpolate_xy_simple(  )
    if (p$spacetime_engine=="LaplacesDemon") pred = spacetime_interpolate_xy_LaplacesDemon(  )
    if (p$spacetime_engine=="inla") pred = spacetime_interpolate_xy_local_inla()

    
   # and then:
  # copy prediction merging/averaging method from xyts  <<<< ----- NOTE TO DO



    if (!is.null(res)) {
      if (exists(p$variogram.engine, res) ) {
        print( "Saving summary statisitics" )
        # save statistics last as this is an indicator of completion of all tasks .. restarts would be broken otherwise
        S[Si,1] = res$varZ
        S[Si,2] = res[[p$variogram.engine]]$varSpatial
        S[Si,3] = res[[p$variogram.engine]]$varObs
        S[Si,4] = res[[p$variogram.engine]]$range
        S[Si,5] = res[[p$variogram.engine]]$phi
        S[Si,6] = res[[p$variogram.engine]]$nu
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

