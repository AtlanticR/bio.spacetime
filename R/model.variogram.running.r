
  model.variogram.running = function( S, v, y, p, empirical.vgm.only=F,
    loc = file.path(  project.datadirectory("snowcrab"), "R", "variograms" ) ) {

    # construct empirical (residual) variogram and model
    yrs =  y + c(-3:3 )
    weights = c( 1/4, 1/2, 3/4, 1, 3/4, 1/2, 1/4 )  # linear decay of weights

    datayrs = unique( S$yr )
    # running variogram 
    out= NULL
    vgm.e = vgm.e0 = buf = np = v.scale = vgm = wts = NULL
    lvars =  c( all.vars( p$kformula) )

    for ( iy in 1:length(yrs) ) { # obtain historical variograms
      if ( yrs[iy] %in% datayrs ) {
        q = which(S$yr==yrs[iy] )
        if (length(q) > 30) {# final check
  #        q = block.xyz2xyz( xyz=q, params=p, engine="R2", variables=lvars )
  #        q$t[!is.finite(q$t)] = mean(q$t, na.rm=T, trim=0.1) # temporaily add  means to force a solution
  #        q$z[!is.finite(q$z)] = mean(q$z, na.rm=T, trim=0.1)
          vgm = try( variogram(id="kv", p$kformula, loc=p$klocs, data=S[q,], boundaries=p$vgm.dist, cutoff=max(p$vgm.dist) + 20, cressie=T ), silent=T)
          if  (! ( "try-error" %in% class(vgm) )  ) { 
            # vgm$gamma = vgm$gamma  # normalised to unit variance for each year
            buf = cbind(buf, vgm$gamma )
            np  = cbind(np, vgm$np)
            wts = c(wts, weights[iy] )
          }
        }
      }
    }
   
    if (is.null(buf)) return(NULL)
    vgm.e = variogram(id="kv", p$kformula, loc=p$klocs, data=S, boundaries=p$vgm.dist, cutoff=max(p$vgm.dist) + 20, cressie=T  ) # obtain template
    vgm.e$gamma = apply( buf, 1, weighted.mean, w=wts, na.rm=T)  # still normalised to unit variance
    vgm.e$np = round( apply( np, 1,  weighted.mean, w=wts, na.rm=T) ) 
    vgm.e$id = "kv"
    
    if (empirical.vgm.only) return (vgm.e)

    vgm.m = try( gstat.model.variogram ( vgm.e, vp=get.variogram.params(v, y) ), silent=T )
    if ( is.null(vgm.m) | ( "try-error" %in% class(vgm.m) ) )  return (NULL )
    

    out = list(vgm.e=vgm.e, vgm.m=vgm.m)
    
    print(vgm.e)
    print(vgm.m)

    if (p$plot.variogram) {
      # plot finalised variograms
      
      fname = file.path(loc, paste(v, y,"png", sep="."))
      dir.create(path=loc, recursive=T, showWarnings=F)
      png(filename=fname)
      plot( x=vgm.e$dist, y=vgm.e$gamma, ylab="Semi-variance", xlab="Distance (km)", pch=20 )
      xaxis = range(p$vgm.dist)
      yaxis = range(vgm.e$gamma)
      textout = paste(
          " Year    = ", y, " ~ mean: ", vyrs.try[vt],
        "\n Kriging = ", deparse(p$kformula),
        "\n Model  = ", vgm.m$model[2],
        "\n Psill     = ", signif(vgm.m$psill[2],3),
        "\n Nugget = ", signif(vgm.m$psill[1],3),
        "\n Range  = ", signif(vgm.m$range[2],3),
        "\n Kappa  = ", signif(vgm.m$kappa[2],3), sep="" )
      lines( variogram.line(vgm.m, maxdist=xaxis[2] ), lty="dashed" )
      text ( x=xaxis[2] *0.5, y= yaxis[2] * 0.6, textout, adj=0, cex=0.8)
      dev.off()
    }

    return( out )
  }


