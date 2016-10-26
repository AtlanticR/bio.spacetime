
spacetime__habitat = function( p, x, pa ) {
   #\\ this is the core engine of spacetime .. localised space-time habiat modelling

  # estimate model parameters
  hmod = try( 
    gam( p$spacetime_engine_modelformula, data=x, weights=Y_wgt, optimizer=c("outer","bfgs")  ) ) 

  if ( "try-error" %in% class(hmod) ) next()
  
  out = try( predict( hmod, newdata=pa, type="response", se.fit=T ) ) 

  if ( "try-error" %in% class( out ) ) return( NULL )

  pa$mean = as.vector(out$fit)
  pa$sd = as.vector(out$se.fit) # this is correct: se.fit== stdev of the mean fit: eg:  https://stat.ethz.ch/pipermail/r-help/2005-July/075856.html

  if (exists( "quantile_bounds", p)) {
    tq = quantile( x[,p$variables$Y], probs=p$quantile_bounds, na.rm=TRUE  )
    bad = which( pa$mean < tq[1] | pa$mean > tq[2]  )
    if (length( bad) > 0) {
      pa$mean[ bad] = NA
      pa$sd[ bad] = NA
    }
  }

  ss = summary(hmod)
  spacetime_stats = list( sdTotal=sd(Y[], na.rm=T), rsquared=ss$r.sq, ndata=ss$n ) # must be same order as p$statsvars

  return( list( predictions=pa, spacetime_stats=spacetime_stats ) )  

}

if (0) {
  Hmodel = NULL
  Hmodel = habitat.model.db( DS="habitat", v=v, yr=y )
  if (is.null( Hmodel)) next()

  Hsim = gam.simulation( M=Hmodel, X= PS, nsims=p$nsims ) #~8min
  rm( Hmodel); gc()

  print("finished H sim")

  oops = which( is.na(Hsim) )
  if (length(oops) > 0)  Hsim[oops ] = 0  # assume to be zero

  Amodel = NULL
  Amodel = habitat.model.db( DS="abundance", v=v, yr=y )
  if (is.null(Amodel)) next()

  Asim = gam.simulation( M=Amodel, X= PS, nsims=p$nsims ) # ~5min
  rm( Amodel); gc()
  print("finished A sim")

  oops = which( is.na(Asim) )
  if (length(oops) > 0)  Asim[oops ] = 0  # assume to be zero

  # Do not extrapolate: trim to XX% quantiles to be a little more conservative
  oopu =  which( Asim > qs[2] )
  if (length(oopu) > 0)  Asim[ oopu ] = qs[2]

  oopl =  which( Asim < qs[1]  )
  if (length(oopl) > 0)  Asim[ oopl ] = 0  # below detection limits

  Asim = Asim * Hsim  # Asim now becomes weighted by Pr of habitat

  Hsim.sa = colSums( Hsim ) # Pr weighted sum of habitat
  totalsurfacearea = mean ( Hsim.sa ) * (p$pres*p$pres)
  totalsurfacearea.sd = sd( Hsim.sa ) * (p$pres*p$pres)
  rm ( Hsim.sa ); gc()

  PS$habitat.mean = apply( Hsim, 1, mean, na.rm=T )
  PS$habitat.sd = apply( Hsim, 1, sd, na.rm=T )

  PS$abundance.mean = apply( Asim, 1, mean, na.rm=T )
  PS$abundance.sd =  apply( Asim, 1, sd, na.rm=T )

  # iAbundance = which ( PS$abundance.mean >= qs[1] )  #  consider abundance only if it is sufficiently precise (low associated variance)
  iHabitat = which( PS$habitat.mean > p$habitat.threshold.quantile  & (PS$habitat.mean - 2 * PS$habitat.sd) > 0 )
  iStations = filter.prediction.locations( DS="limit.to.near.survey.stations", PS=PS, y=y, p=p )  # consider locations only if close to a survey location (minimal extrapolation)
  # iHA = intersect(iHabitat, iAbundance)  # good locations for habitat and abundance prediction
  iE = union( iStations, intersect(iHabitat, iStations ))  #  good locations for habitat and abundance prediction, but filtered for proximity to survey stations
  rm( iStations); gc()

  # levelplot( habitat.mean ~ plon + plat, data=PS[,], aspect="iso" )
  # levelplot( habitat.mean ~ plon + plat, data=PS[iE,], aspect="iso" )

  fn.res = file.path( loc.sol, paste( "PS.simulation.means", v, y, "rdata", sep="." ) )
  print (fn.res )
  save( PS, file=fn.res, compress=T )

  K = NULL

  fn.K = file.path( loc.res, paste( "K", v, y, "rdata", sep="." ) )
  print(fn.K)

  PS$abundance.mean.log = log10( PS$abundance.mean )  # only used for plotting
  er = log10( qs  )

  PS$abundance.mean.log [ which( PS$abundance.mean.log < er[1] ) ] = er[1]
  # levelplot( abundance.mean.log ~ plon + plat, data=PS, aspect="iso" )

  if ("map.habitat" %in% p$ItemsToMap ) {
    datarange = seq( 0, 1, length.out=150)
    cols = color.code( "seis", datarange )
    outfn = paste( "prediction.habitat.mean", v, y, sep=".")
    map( xyz=PS[,c("plon", "plat", "habitat.mean")], cfa.regions=T, depthcontours=T, pts=NULL, annot=paste( v, y ),
      annot.cex=p$annot.cex, corners=p$planar.corners, fn=outfn, loc=loc.map, at=datarange,
      col.regions=cols, rez=c(p$pres,p$pres) )
  }

  if ("map.abundance" %in% p$ItemsToMap ) {
    datarange = seq( er[1], er[2], length.out=150)
    cols = color.code( "seis", datarange )
    outfn = paste( "prediction.abundance.mean", v, y, sep=".")
    map( xyz=PS[ , c("plon", "plat", "abundance.mean.log")], cfa.regions=T, depthcontours=T, pts=NULL, annot= paste( v, y ),
      annot.cex=p$annot.cex, corners=p$planar.corners, fn=outfn, loc=loc.map, at=datarange,
      col.regions=cols, rez=c(p$pres,p$pres) )
  }

  if ("map.abundance.estimation" %in% p$ItemsToMap ) {
    tomap = NULL
    tomap = PS[ , c("plon", "plat", "abundance.mean.log")]
    tomap$abundance.mean.log[ setdiff( 1:nrow(PS), iE ) ] = log10( er[1] ) - 1  # just a value smaller than the lower bound
    datarange = seq( er[1], er[2], length.out=150)
    cols = color.code( "seis", datarange )
    outfn = paste( "prediction.abundance.mean.estimationarea", v, y, sep=".")
    map( xyz=tomap, cfa.regions=T, depthcontours=T, pts=NULL,
      annot=paste( v, y ), annot.cex=p$annot.cex, corners=p$planar.corners, fn=outfn, loc=loc.map, at=datarange,
      col.regions=cols, rez=c(p$pres,p$pres) )
    rm (tomap) ; gc()
  }

  PS$abundance.mean.log = NULL
  gc()

  print( "finished maps")

  for (r in p$regions ){

    iRegion.PS = filter.region.polygon(x=PS[ , c("plon", "plat")], region=r, planar=T)
    iEastof250 = which( PS$plon > 250 )
    iRegion.PS = intersect( iRegion.PS, iEastof250 )

    iHabitatSubarea = intersect( iRegion.PS, iHabitat ) # plotting surface and real habitat area
    iEstimationArea = intersect( iRegion.PS, iE )

    if ( length( iEstimationArea ) > 10 ) {

      hhsum = colSums(  Hsim[ iHabitatSubarea , ] ) # area weighted by Pr
      sa.region = mean( hhsum ) * (p$pres*p$pres)
      sa.sd = sd( hhsum ) * (p$pres*p$pres)

      hhsum = colSums( Hsim[ iEstimationArea , ] )
      sa.estim = mean( hhsum ) * (p$pres*p$pres)
      sa.estim.sd = sd( hhsum ) * (p$pres*p$pres)

      aa.sums = apply( Asim[ iEstimationArea,  ] , 2, sum ) # abundance weighted by Pr
      V.mean = mean( aa.sums )
      V.sd = sd( aa.sums )
      V.sd.ln = sd( log( aa.sums ), na.rm=T )
      ci = quantile( aa.sums, probs=c(0.025, 0.5, 0.975), na.rm=T,names=F )

    } else {
      V.mean = NA
      V.sd = NA
      V.sd.ln = NA
      sa.estim = NA
      sa.estim.sd = NA
      sa.region = NA
      sa.sd = NA
      ci = rep( NA, 3 )
    }

    bb.sums = apply( Asim[ iHabitatSubarea,  ] , 2, sum ) # abundance weighted by Pr
    W.mean = mean( bb.sums )
    W.sd = sd( bb.sums )
    W.ci = quantile( bb.sums, probs=c(0.025, 0.5, 0.975), na.rm=T,names=F )

    L = data.frame( yr=y, vars=v, region=r,
        total=V.mean, total.sd=V.sd, total.sd.ln=V.sd.ln, median=ci[2], lbound=ci[1], ubound=ci[3],
        ss=W.mean, ss.sd=W.sd, ss.median=W.ci[2], ss.lbound=W.ci[1], ss.ubound=W.ci[3],
        sa.total=totalsurfacearea, sa.total.sd=totalsurfacearea.sd,
        sa.region=sa.region, sa.region.sd=sa.sd,
        sa.estimation=sa.estim, sa.estimation.sd=sa.estim.sd,

        datestamp=as.character( Sys.time() ) )
    K = rbind( K, L )
  } # regions

  save(K, file=fn.K, compress=T)
}
