
spacetime_prediction_area = function(p, Si, dist.cur ){
  # So, YiU and dist.cur determine the data entering into local model construction
  # but predictions are expensive, esp in space-time and so a smaller spatial raange is selected
  # this is determined by p$spacetime.prediction.dist.min
  dist_prediction = min( dist.cur*p$spacetime.prediction.range.proportion,  p$spacetime.prediction.dist.min )
  windowsize.half = floor(dist_prediction/p$pres) # convert distance to discretized increments of row/col indices
  pa_w = -windowsize.half : windowsize.half
  pa_w_n = length(pa_w)
  iwplon = p$rcS[Si,1] + pa_w
  iwplat = p$rcS[Si,2] + pa_w

  pa = data.frame( iplon = rep.int(iwplon, pa_w_n) , 
                   iplat = rep.int(iwplat, rep.int(pa_w_n, pa_w_n)) )

  bad = which( (pa$iplon < 1 & pa$iplon > p$nplons) | (pa$iplat < 1 & pa$iplat > p$nplats) )
  if (length(bad) > 0 ) pa = pa[-bad,]
  if (nrow(pa)< 5) return(NULL)

  pc_rc = paste( pa$iplon, pa$iplat, sep="~" )
  pa$i = match( pc_rc, p$rcP$rc)
  bad = which( !is.finite(pa$i))
  if (length(bad) > 0 ) pa = pa[-bad,]

  pa_n = nrow(pa)
  if ( pa_n < 5) return(NULL)

  Ploc = switch( p$storage.backend, 
    bigmemory.ram=attach.big.matrix(p$ptr$Ploc), 
    bigmemory.filebacked=attach.big.matrix(p$ptr$Ploc), 
    ff=p$ptr$Ploc )

  pa$plon = Ploc[ pa$i, 1]
  pa$plat = Ploc[ pa$i, 2]

  # prediction covariates i.e., independent variables/ covariates
  if (exists("COV", p$variables)) {

    Pcov = switch( p$storage.backend, 
      bigmemory.ram=attach.big.matrix(p$ptr$Pcov), 
      bigmemory.filebacked=attach.big.matrix(p$ptr$Pcov), 
      ff=p$ptr$Pcov )

    for (ci in 1:length(p$variables$COV)) {
      pa[,p$variables$COV[ci]] = Pcov[ pa$i, ci ]
    }
  }

  if (0) {
    Sloc = switch( p$storage.backend, 
      bigmemory.ram=attach.big.matrix(p$ptr$Sloc), 
      bigmemory.filebacked=attach.big.matrix(p$ptr$Sloc), 
      ff=p$ptr$Sloc )
    Yloc = switch( p$storage.backend, 
      bigmemory.ram=attach.big.matrix(p$ptr$Yloc), 
      bigmemory.filebacked=attach.big.matrix(p$ptr$Yloc), 
      ff=p$ptr$Yloc )

    plot( Yloc[U,1]~ Yloc[U,2], col="red", pch=".") # all data
    points( Yloc[YiU,1] ~ Yloc[YiU,2], col="green" )  # with covars and no other data issues
    points( Sloc[Si,1] ~ Sloc[Si,2], col="blue" ) # statistical locations
    points( p$plons[p$rcS[Si,1]] ~ p$plats[p$rcS[Si,2]] , col="purple", pch=25, cex=2 ) # check on p$rcS indexing
    points( p$plons[pa$iplon] ~ p$plats[ pa$iplat] , col="cyan", pch=".", cex=0.01 ) # check on Proc iplat indexing
    points( Ploc[pa$i,1] ~ Ploc[ pa$i, 2] , col="black", pch=20, cex=0.7 ) # check on pa$i indexing -- prediction locations
  }
  
  if ( exists("COV", p$variables) ) {
    pvars = c("plon", "plat", "i", p$variables$COV)
  } else {
    pvars = c("plon", "plat", "i")
  }

  if ( ! exists("TIME", p$variables) ) {
    pa = pa[, pvars]
  } else {
    Ptime = switch( p$storage.backend, 
      bigmemory.ram=attach.big.matrix(p$ptr$Ptime), 
      bigmemory.filebacked=attach.big.matrix(p$ptr$Ptime), 
      ff=p$ptr$Ptime )
    pa = cbind( pa[ rep.int(1:pa_n, length(Ptime)), pvars ], 
                     rep.int(Ptime[], rep(pa_n,length(Ptime) )) )
    names(pa) = c(pvars, "tiyr" )
  }
  return(pa)

}


