

spacetime = function( p, DATA, OUT=NULL, overwrite=NULL, DS=NULL, method="inla" ) {
	#\\ localized modelling of space and time data to predict/interpolate upon a grid OUT

	if (!is.null(DS)) {
	#\\ data extraction layer

      if ( DS=="spatial.covariance")  {
        stats = NULL
        if (file.exists( p$fn.results.covar) ) load( p$fn.results.covar )
        return(stats)
      }

      if ( DS =="inla.predictions" ) {
        preds = NULL
        if (file.exists( p$fn.P ) ) load( p$fn.P )
        return( preds )
      }

      if ( DS =="inla.statistics" ) {
        stats = NULL
        if (file.exists( p$fn.S) ) load( p$fn.S )
        return( stats )
      }

	}


    # set up the data and problem using bigmemory data objects
    print( paste( "Temporary files are being created at:", p$tmp.datadir ) )
    p = spacetime.db( p=p, DS="bigmemory.filenames" )

    dir.create( p$rootdir, showWarnings=FALSE, recursive =TRUE)

    if (is.null(overwrite)) {
	    tmpfiles = spacetime.db( p=p, DS="bigmemory.filelist" )
		for (tf in tmpfiles) {
	  		if (file.exists( tf)) {
	  			cat( "Temporary files exist from a previous run found. \n")
	  			cat( "Send explicit overwrite=TRUE or overwrite=FALSE to proceed. \n") 
	  			stop()
	  		}
	  	}
	}

	if (method=="spatial.covariance" ) {

		if (is.null(overwrite) || overwrite) {
	      # init input data objects
	      spacetime.db( p=p, DS="bigmemory.inputs.data", B=DATA )
	      # init output data objects
	      spacetime.db( p=p, DS="statistics.bigmemory.initialize" )
		}

  	  # define boundary polygon for data .. 
  	  # this trims the prediction/statistics locations to speed things up a little ..
	  spacetime.db( p, DS="boundary.redo" ) # ~ 5 min

      if (0) {
        # to reset results manually .. just a template
        # p = spacetime.db( p=p, DS="bigmemory.filenames" )
        S = bigmemory::attach.big.matrix(p$descriptorfile.S, path=p$tmp.datadir)  # statistical outputs
        hist(S[,1] )
        o = which( S[,1] > xxx )
        S[o,] = NA
        S[sS$problematic,] = NA
        o = which( S[,1] < yyy )
        S[o,] = NA
      }

      sS = spacetime.db( p, DS="statistics.bigmemory.status" )
      sS$n.incomplete / ( sS$n.problematic + sS$n.incomplete + sS$n.complete)

      p = make.list( list( jj=sample( sS$incomplete ) ), Y=p ) # random order helps use all cpus
      parallel.run( spacetime.covariance.spatial, p=p ) # no more GMT dependency! :)
      # spacetime.covariance.spatial( p=p )  # if testing serial process

      # save to file
      print( paste( "Results are being saved to:", p$fn.results.covar ) )
      stats = bigmemory::attach.big.matrix(p$descriptorfile.S, path=p$tmp.datadir)  # statistical outputs
      stats = as.data.frame( stats[] )
      save(stats, file=p$fn.results.covar, compress=TRUE )

      print( paste( "Temporary files are being deleted at:", p$tmp.datadir, "tmp" ) )
      spacetime.db( p=p, DS="bigmemory.cleanup" )

      return( p )
    }


    # ------------------------

    
    if (  method =="inla.interpolations" ) {

		if (is.null(overwrite) || overwrite) {
	      # init input data objects
	      spacetime.db( p=p, DS="bigmemory.inputs.data", B=DATA )
	      spacetime.db( p=p, DS="bigmemory.inputs.prediction", B=OUT) # covas on prediction locations
	      # init output data objects
	      spacetime.db( p=p, DS="statistics.bigmemory.initialize" )
	      spacetime.db( p=p, DS="predictions.bigmemory.initialize" )
		}

      # run the beast .. warning this will take a very long time! (weeks)
      sS = spacetime.db( p, DS="statistics.bigmemory.status" )
      sS$n.incomplete / ( sS$n.problematic + sS$n.incomplete + sS$n.complete)

      p = make.list( list( jj=sample( sS$incomplete ) ), Y=p ) # random order helps use all cpus
      parallel.run( spacetime.interpolate.inla.local, p=p ) # no more GMT dependency! :)
      # spacetime.interpolate.inla.local( p=p, debugrun=TRUE )  # if testing serial process

      if (0) {
        # for checking status of outputs **during** parallel runs: they access bigmemory temporary files
#        bathymetry.figures( DS="statistics", p=p )
#        bathymetry.figures( DS="predictions", p=p )
#        bathymetry.figures( DS="predictions.error", p=p )
        p = spacetime.db( p=p, DS="bigmemory.filenames" )
        S = bigmemory::attach.big.matrix(p$descriptorfile.S, path=p$tmp.datadir)  # statistical outputs
        hist(S[,1] )
        o = which( S[,1] > 600 )
        S[o,] = NA
        S[sS$problematic,] = NA
        o = which( S[,1] < 10 )
        S[o,] = NA
      }

      # save to file
    
        pp = bigmemory::attach.big.matrix(p$descriptorfile.P, path=p$tmp.datadir)  # predictions
        preds = pp[]
        ppl = bigmemory::attach.big.matrix(p$descriptorfile.Ploc, path=p$tmp.datadir)
        predloc = ppl[]
        preds = as.data.frame( cbind ( predloc, preds ) )
        names(preds) = c( "plon", "plat", "ndata", "mean", "sdev" )
        save( preds, file=p$fn.P, compress=TRUE )


	  # this also rescales results to the full domain
        #\\ statistics are stored at a different resolution than the final grid
        #\\   this fast interpolates the solutions to the final grid
        S = bigmemory::attach.big.matrix(p$descriptorfile.S, path=p$tmp.datadir)  # statistical outputs
        ss = as.data.frame( S[] )
        statnames0 = c( "range", "range.sd", "spatial.var", "observation.var"  )
        statnames  = c( "range", "range.sd", "spatial.sd", "observation.sd"  )
        datalink   = c( "log", "log", "log", "log" )  # a log-link seems appropriate for these data
        names(ss) = statnames0
        ssl = bigmemory::attach.big.matrix(p$descriptorfile.Sloc, path=p$tmp.datadir)  # statistical output locations
        sslocs = as.data.frame(ssl[]) # copy
        names(sslocs) = p$variables$LOCS
        ss = cbind( sslocs, ss )
        rm (S)
        ss$spatial.sd = sqrt( ss$spatial.var )
        ss$observation.sd = sqrt( ss$observation.var )
        ss$spatial.var = NULL
        ss$observation.var = NULL

        # trim quaniles in case of extreme values
        for ( v in statnames ) {
          vq = quantile( ss[,v], probs= c(0.025, 0.975), na.rm=TRUE )
          ii = which( ss[,v] < vq[1] )
          if ( length(ii)>0) ss[ii,v] = vq[1]
          jj = which( ss[,v] > vq[2] )
          if ( length(jj)>0) ss[jj,v] = vq[2]
        }

        locsout = expand.grid( p$plons, p$plats ) # final output grid
        attr( locsout , "out.attrs") = NULL
        names( locsout ) = p$variables$LOCS

        stats = matrix( NA, ncol=length(statnames), nrow=nrow( locsout) )  # output data
        colnames(stats)=statnames

        for ( iv in 1:length(statnames) ) {
          vn = statnames[iv]
          # create a "surface" and interpolate to larger grid using
          # (gaussian) kernel-based smooth on the log-scale
          z = log( matrix( ss[,vn], nrow=length(p$sbbox$plons), ncol=length( p$sbbox$plats) ) )
          RES = NULL

          interp.method = "kernel.density" # default 

          if (interp.method == "kernel.density") {
	          # create a "surface" and interpolate to a different grid using (gaussian) kernel-based smooth via FFT
		    require(fields)
		  	isurf = fields::interp.surface( list( x=p$sbbox$plons, y=p$sbbox$plats, z=z), loc=locsout  )
		  	zout = matrix( isurf, nrow=length(p$plons), ncol=length( p$plats) )
		  	RES = fields::image.smooth( zout, theta=p$dist.mwin, xwidth=p$dist.mwin*10, ywidth=p$dist.mwin*10 ) # 10 SD of the normal kernel
            if ( !is.null( RES )) stats[,iv] = exp( RES$z ) # return to correct scale
          }

          if (interp.method=="inla.fast") { # fast, but not fast enough for prime time yet
            # interpolation using inla is also an option
            # but will require a little more tweaking as it was a bit slow
            range0 = median( ss$range, na.rm=TRUE )
            oo = which( is.finite( ss[,vn] ) )
            if ( length(oo) < 30 ) next()
            RES = spacetime.interpolate.inla.singlepass (
              Y=ss[oo,vn], locs=ss[oo, p$variables$LOCS], plocs=locsout, method="fast", link=datalink[iv] )
            if ( !is.null( RES )) stats[,iv] = RES$xmean
          }
        }

        save( stats,  file=p$fn.S, compress=TRUE )


        plotdata=FALSE ## to debug
        if (plotdata) {
          p$spatial.domain="canada.east"  # force isobaths to work in levelplot
          datarange = log( c( 5, 1200 ))
          dr = seq( datarange[1], datarange[2], length.out=150)
          oc = landmask( db="worldHires", regions=c("Canada", "US"),
                         return.value="not.land", tag="predictions" )  ## resolution of "predictions" which is the final grid size
          toplot = cbind( locsout, z=(stats[,"range"]) )[oc,]
          resol = c(p$dist.mwin,p$dist.mwin)
          levelplot( log(z) ~ plon + plat, toplot, aspect="iso", at=dr, col.regions=color.code( "seis", dr) ,
            contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE), cex=2, resol=resol,
            panel = function(x, y, subscripts, ...) {
              panel.levelplot (x, y, subscripts, aspect="iso", rez=resol, ...)
              cl = landmask( return.value="coast.lonlat",  ylim=c(36,53), xlim=c(-72,-45) )
              cl = lonlat2planar( data.frame( cbind(lon=cl$x, lat=cl$y)), proj.type=p$internal.crs )
              panel.xyplot( cl$plon, cl$plat, col = "black", type="l", lwd=0.8 )
            }
          )
          p$spatial.domain="canada.east.highres"
        }
      }

      # clean up bigmemory files
      spacetime.db( p=p, DS="bigmemory.cleanup" )
      return(p)
  }


}



