
spacetime_interpolate = function( ip=NULL, p ) {
  #\\ core function to intepolate (model and predict) in parllel

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  #---------------------
  # data for modelling
  # dependent vars # already link-transformed in spacetime_db("dependent")
  S = spacetime_attach( p$storage.backend, p$ptr$S )
  Sflag = spacetime_attach( p$storage.backend, p$ptr$Sflag )
  
  # force copy into RAM
  Sloc = spacetime_attach( p$storage.backend, p$ptr$Sloc )[]
  Yloc = spacetime_attach( p$storage.backend, p$ptr$Yloc )[]

  Yi = spacetime_attach( p$storage.backend, p$ptr$Yi )
  Yi = as.vector(Yi[])  #force copy to RAM

  # main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    Si = p$runs[ iip, "locs" ]
    if ( is.infinite( Sflag[Si] ) ) next() 
    if ( !is.nan( Sflag[Si] ) ) next() 
    Sflag[Si] = Inf   # over-written below if successful else if a run fails it does not get revisited 
    print( iip )

    # find data withing a given distance / number 
    pib = point_in_block( Sloc=Sloc, Si=Si, Yloc=Yloc, Yi=Yi, 
      dist.max=p$dist.max, dist.min=p$dist.min, n.min=p$n.min, n.max=p$n.max, 
      upsampling=p$upsampling, downsampling=p$downsampling, resize=TRUE ) 
    
    if ( is.null(pib)) next()
    
    dist.cur = pib$dist
    ndata = length(pib$U)
    if ((ndata < p$n.min) | (ndata > p$n.max) ) next()
    YiU = Yi[pib$U]  
    rm(pib)
    
    # construct prediction/output grid area ('pa')
    pa = NULL
    pa = spacetime_prediction_area( p, Si, dist.cur ) 
    if (is.null(pa)) next()
    
    res = NULL
    res = spacetime_model_predict( p, Si, YiU, pa )     # model and prediction
    if ( is.null(res)) next()
    rm(pa)
    
    for ( k in 1: length(p$statsvars) ) {
      if (exists( p$statsvars[k], res$spacetime_stats )) {
        S[Si,k] = res$spacetime_stats[[ p$statsvars[k] ]]
      }
    }

    spacetime_predictions_save( p, res$predictions ) # update P (predictions) .. slow!! .. try diff cache size for P, Pn Psd

      if (0) {
        v = res$predictions
        if ( exists("TIME", p$variables) ){
          v = v[which( v[,p$variables$TIME]==1990.55),]
        }
        require(lattice)
        levelplot( mean ~ plon+plat, v, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
      }

   
    # ----------------------
    # save statistics: do last. it is an indicator of completion of all tasks 
    # .. restarts would be broken otherwise
    Sflag[Si] = 1  # done .. any finite value

  }  # end for loop

}



