
spacetime_interpolate = function( ip=NULL, p ) {

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  #---------------------
  # data for modelling
  # dependent vars # already link-transformed in spacetime_db("dependent")
  Yi = p$ptr$Yi
  Y =  p$ptr$Y 
  Yloc = p$ptr$Yloc 
  Sloc = p$ptr$Sloc 

  # main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    Si = p$runs[ iip, "locs" ]
    S = p$ptr$S  # statistical outputs inside loop to safely save data and pass onto other processes
    if ( is.infinite( S[Si,1] ) ) next() 
    if ( !is.nan( S[Si,1] ) ) next() 
    # Si = 31133  problem
    S[Si,1] = Inf   # over-written below if successful else if a run fails it does not get revisited 
    close(S)
    print( iip )

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
    
    # construct prediction/output grid area ('pa')
    pa = spacetime_prediction_area( p, Si, dist_cur ) 
    res = spacetime_model_predict( p, Si, YiU, pa )     # model and prediction
    if ( is.null(res)) next()
    spacetime_predictions_save( p, res$predictions ) # update P (predictions) .. slow!! .. try diff cache size for P, Pn Psd

      if (0) {
        v = res$predictions
        if ( exists("TIME", p$variables) ){
          v = v[which( v[,p$variables$TIME]==1990.55),]
        }
        require(lattice)
        levelplot( mean ~ plon+plat, v, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
      }
    rm(res, pa); gc()

    # ----------------------
    # save statistics: do last. it is an indicator of completion of all tasks 
    # .. restarts would be broken otherwise
    S = ( p$ptr$S )
    for ( k in 1: length(p$statsvars) ) {
      if (exists( p$statsvars[k], res$spacetime_stats )) {
        S[Si,k] = res$spacetime_stats[[ p$statsvars[k] ]]
      }
    }
    close(S)

  }  # end for loop

}



