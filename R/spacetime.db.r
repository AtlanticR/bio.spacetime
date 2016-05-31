
  spacetime.db = function( DS, p, B=NULL, grp=NULL ) {
    #// usage: low level function to convert data into bigmemory obects to permit parallel
    #// data access and manipulation
    #// B is the xyz  or xytz data to work upon

    if (DS %in% "bigmemory.filenames" ) {

      # create file backed bigmemory objects

      p$tmp.datadir = file.path( p$project.root, "tmp" )
      if( !file.exists(p$tmp.datadir)) dir.create( p$tmp.datadir, recursive=TRUE, showWarnings=FALSE )

      # input data stored as a bigmatrix to permit operations with min memory usage
      # split into separate components to minimize filelocking conflicts
      p$backingfile.Y = "input.Y.bigmatrix.tmp"
      p$descriptorfile.Y = "input.Y.bigmatrix.desc"

      p$backingfile.X = "input.X.bigmatrix.tmp"
      p$descriptorfile.X = "input.X.bigmatrix.desc"

      p$backingfile.LOCS = "input.LOCS.bigmatrix.tmp"
      p$descriptorfile.LOCS = "input.LOCS.bigmatrix.desc"

      p$backingfile.P = "predictions.bigmatrix.tmp"
      p$descriptorfile.P = "predictions.bigmatrix.desc"

      p$backingfile.S = "statistics.bigmatrix.tmp"
      p$descriptorfile.S = "statistics.bigmatrix.desc"

      p$backingfile.Pcov = "predictions_cov.bigmatrix.tmp"
      p$descriptorfile.Pcov = "predictions_cov.bigmatrix.desc"

      p$backingfile.Ploc = "predictions_loc.bigmatrix.tmp"
      p$descriptorfile.Ploc = "predictions_loc.bigmatrix.desc"

      p$backingfile.Sloc = "statistics_loc.bigmatrix.tmp"
      p$descriptorfile.Sloc = "statistics_loc.bigmatrix.desc"


      return(p)
    }

    # --------------------------

    if (DS %in% "bigmemory.cleanup" ) {
      # load bigmemory data objects pointers
      p = spacetime.db( p=p, DS="bigmemory.filenames" )
      todelete = file.path( p$tmp.datadir,
        c( p$backingfile.P, p$descriptorfile.P,
           p$backingfile.S, p$descriptorfile.S,
           p$backingfile.Sloc, p$descriptorfile.Sloc,
           p$backingfile.Ploc, p$descriptorfile.Ploc,
           p$backingfile.Pcov, p$descriptorfile.Pcov,
           p$backingfile.Y, p$descriptorfile.Y,
           p$backingfile.X, p$descriptorfile.X,
           p$backingfile.LOCS, p$descriptorfile.LOCS
      ))
      for (fn in todelete ) if (file.exists(fn)) file.remove(fn)
      return( todelete )
    }

    # ------------------

    if (DS == "bigmemory.inputs.data" ) {
      spacetime.db( p=p, DS="bigmemory", B=B, grp="dependent" )
      spacetime.db( p=p, DS="bigmemory", B=B, grp="coordinates" )
      spacetime.db( p=p, DS="bigmemory", B=B, grp="covariates" )
    }

    # ------------------

    if (DS == "bigmemory.inputs.prediction" ) {
      spacetime.db( p=p, DS="bigmemory", B=B, grp="prediction.coordinates" )
      spacetime.db( p=p, DS="bigmemory", B=B, grp="prediction.covariates" )
    }


    # ------------------

    if (DS == "bigmemory" ) {
      # create file backed bigmemory objects
      p = spacetime.db( p=p, DS="bigmemory.filenames" )  # load bigmemory data objects pointers

      if (grp=="dependent") {
        # dependent variable
        fn.Y = file.path(p$tmp.datadir, p$backingfile.Y )
        if ( file.exists( fn.Y) ) file.remove( fn.Y)
        Y = filebacked.big.matrix( nrow= nrow(B), ncol=1, type="double", dimnames=NULL, separated=FALSE,
          backingpath=p$tmp.datadir, backingfile=p$backingfile.Y, descriptorfile=p$descriptorfile.Y )
        if ( "data.frame" %in% class(B) ) {
          Y[] = as.matrix( B[ , p$variables$Y ] )
        } else if ( "SpatialGridDataFrame" %in% class(B) ) {
          Y[] = as.matrix( slot(B, "data")[, p$variables$Y ]  )
        }
        if ( exists("spacetime.link", p)) Y[] = p$spacetime.link ( Y[] )
      }

      if (grp=="covariates") {
        # independent variables/ covariates
        if ( exists( "X", p$variables) ) {
          fn.X = file.path(p$tmp.datadir, p$backingfile.X )
          if ( file.exists( fn.X) ) file.remove( fn.X)
          X = filebacked.big.matrix( nrow=nrow(B), ncol=length( p$variables$X ), type="double", dimnames=NULL, separated=FALSE,
            backingpath=p$tmp.datadir, backingfile=p$backingfile.X, descriptorfile=p$descriptorfile.X )
          if ( "data.frame" %in% class(B) ) {
            X[] = as.matrix( B[ , p$variables$X ] )
          } else if ( "SpatialGridDataFrame" %in% class(B) ) {
            X[] = as.matrix( slot(B, "data")[, p$variables$X ]  )
          }
        }
      }

      if (grp=="coordinates") {
        # coordinates
        fn.LOC = file.path(p$tmp.datadir, p$backingfile.LOC )
        if ( file.exists( fn.LOC) ) file.remove( fn.LOC)
        LOCS = filebacked.big.matrix( nrow=nrow(B), ncol=2, type="double", dimnames=NULL, separated=FALSE,
            backingpath=p$tmp.datadir, backingfile=p$backingfile.LOCS, descriptorfile=p$descriptorfile.LOCS )
        if ( "data.frame" %in% class(B) ) {
          LOCS[] = as.matrix( B[ , p$variables$LOCS ] )
        } else if ( "SpatialGridDataFrame" %in% class(B) ) {
          LOCS[] = as.matrix( coordinates(B) )
        }
      }

      if (grp=="prediction.coordinates") {
        # prediction coordinates
        fn.Ploc = file.path(p$tmp.datadir, p$backingfile.Ploc )
        if ( file.exists( fn.Ploc) ) file.remove( fn.Ploc )
        Ploc = filebacked.big.matrix( nrow=nrow(B), ncol=2, type="double", dimnames=NULL, separated=FALSE,
           backingpath=p$tmp.datadir, backingfile=p$backingfile.Ploc, descriptorfile=p$descriptorfile.Ploc )
        if ( "data.frame" %in% class(B) ) {
          Ploc[] = as.matrix( B[ , p$variables$LOCS ] )
        } else if ( "SpatialGridDataFrame" %in% class(B) ) {
          Ploc[] = as.matrix( coordinates(B) )
        }
      }

      if (grp=="prediction.covariates") {
        # prediction covariates i.e., independent variables/ covariates
        if ( exists( "X", p$variables) ) {
          fn.Pcov = file.path(p$tmp.datadir, p$backingfile.Pcov )
          if ( file.exists( fn.Pcov) ) file.remove( fn.Pcov)
          Pcov = filebacked.big.matrix( nrow=nrow(B), ncol=length( p$variables$X ), type="double", dimnames=NULL, separated=FALSE,
            backingpath=p$tmp.datadir, backingfile=p$backingfile.Pcov, descriptorfile=p$descriptorfile.Pcov )
          if ( "data.frame" %in% class(B) ) {
            Pcov[] = as.matrix( B[ , p$variables$X ] )
          } else if ( "SpatialGridDataFrame" %in% class(B) ) {
            Pcov[] = as.matrix( slot(B, "data")[, p$variables$X ]  )
          }
        }
      }

      if (grp=="statistics.coordinates") {
        # statistics coordinates
        fn.Sloc = file.path(p$tmp.datadir, p$backingfile.Sloc )
        if ( file.exists( fn.Sloc) ) file.remove( fn.Sloc )
        coords = expand.grid( p$sbbox$plons, p$sbbox$plats )
        Sloc = filebacked.big.matrix( nrow=nrow(coords), ncol=2, type="double", dimnames=NULL, separated=FALSE,
           backingpath=p$tmp.datadir, backingfile=p$backingfile.Sloc, descriptorfile=p$descriptorfile.Sloc )
        Sloc[] = as.matrix( coords )
      }

      if (grp=="statistics.results") {
        # statistics results output file .. initialize
        coords = expand.grid( p$sbbox$plons, p$sbbox$plats )

        fn.S = file.path(p$tmp.datadir, p$backingfile.S )
        if ( file.exists( fn.S) ) file.remove( fn.S)
        S = filebacked.big.matrix( nrow=nrow(coords), ncol= length( p$statsvars ), type="double", init=NA, dimnames=NULL, separated=FALSE,
          backingpath=p$tmp.datadir, backingfile=p$backingfile.S, descriptorfile=p$descriptorfile.S )
      }

      return( "complete" )
    }



    # ----------------
    if (DS %in% c( "predictions", "predictions.redo", "predictions.bigmemory.initialize" )  ) {
      # load bigmemory data objects pointers for predictions
      p = spacetime.db( p=p, DS="bigmemory.filenames" )
      rootdir = file.path( p$project.root, "interpolated" )
      dir.create( rootdir, showWarnings=FALSE, recursive =TRUE)
      fn.P =  file.path( rootdir, paste( "spacetime", "predictions", p$spatial.domain, "rdata", sep=".") )
      if ( DS=="predictions" ) {
        preds = NULL
        if (file.exists( fn.P ) ) load( fn.P )
        return( preds )
      }
      if ( DS=="predictions.bigmemory.initialize" ) {
        # predictions storage matrix (discretized)
        fn.P = file.path(p$tmp.datadir, p$backingfile.P )
        if ( file.exists( fn.P) ) file.remove( fn.P)
        # contains c(count, pred.mean, pred.sd)
        P = filebacked.big.matrix( nrow=p$nplon * p$nplat, ncol=3, type="double", init=NA, dimnames=NULL, separated=FALSE,
          backingpath=p$tmp.datadir, backingfile=p$backingfile.P, descriptorfile=p$descriptorfile.P )
        return( fn.P )
      }

      if ( DS =="predictions.redo" ) {
        pp = bigmemory::attach.big.matrix(p$descriptorfile.P, path=p$tmp.datadir)  # predictions
        preds = pp[]
        ppl = bigmemory::attach.big.matrix(p$descriptorfile.Ploc, path=p$tmp.datadir)
        predloc = ppl[]
        preds = as.data.frame( cbind ( predloc, preds ) )
        names(preds) = c( "plon", "plat", "ndata", "mean", "sdev" )
        save( preds, file=fn.P, compress=TRUE )
        return(fn.P)
      }
    }

    # -----------------

    if (DS == "statistics.box")  {
      sbbox = list( plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$dist.mwin ),
                    plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$dist.mwin )
      )
      return(sbbox)
    }

    # -----------------

    if (DS %in% c( "boundary.redo", "boundary" ) )  {

      p = spacetime.db( p=p, DS="bigmemory.filenames" )
      fn =  file.path(p$tmp.datadir, "boundary.rdata" )

      if (DS=="boundary") {
        boundary = NULL
        if( file.exists(fn)) load( fn)
        return( boundary )
      }

      # load bigmemory data objects pointers

      # data:
      Y = bigmemory::attach.big.matrix(p$descriptorfile.Y, path=p$tmp.datadir )
      LOCS = bigmemory::attach.big.matrix(p$descriptorfile.LOCS, path=p$tmp.datadir )
      hasdata = 1:length(Y)
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA
      # covariates (independent vars)
      if ( exists( "X", p$variables) ) {
        X = bigmemory::attach.big.matrix(p$descriptorfile.X, path=p$tmp.datadir )
        if ( length( p$variables$X ) == 1 ) {
          bad = which( !is.finite( X[]) )
        } else {
          bad = which( !is.finite( rowSums(X[])) )
        }
        if (length(bad)> 0 ) hasdata[bad] = NA
      }
      ii = na.omit(hasdata)
      ndata = length(ii)
      locs_noise = LOCS[ii,] + runif( ndata*2, min=-p$pres*p$spacetime.noise, max=p$pres*p$spacetime.noise )
      maxdist = max( diff( range( LOCS[ii,1] )), diff( range( LOCS[ii,2] )) )

      convex = -0.04
      if (exists( "mesh.boundary.convex", p) ) convex=p$mesh.boundary.convex
      resolution = 125
      if (exists( "mesh.boundary.resolution", p) ) resolution=p$mesh.boundary.resolution
      boundary=list( polygon = inla.nonconvex.hull(  LOCS[ii,], convex=convex, resolution=resolution ) )

      Sloc = bigmemory::attach.big.matrix(p$descriptorfile.Sloc , path=p$tmp.datadir )  # statistical output locations
      boundary$inside.polygon = point.in.polygon( Sloc[,1], Sloc[,2],
          boundary$polygon$loc[,1], boundary$polygon$loc[,2], mode.checked=TRUE)

      save( boundary, file=fn, compress=TRUE )

      plot( LOCS[ii,], pch="." ) # data locations
      lines( boundary$polygon$loc , col="green" )


      return( fn )
    }


    # -----------------

    if (DS %in% c( "statistics", "statistics.redo", "statistics.bigmemory.initialize",
                   "statistics.bigmemory.size" , "statistics.bigmemory.status"  )  ) {

      # load bigmemory data objects pointers
      p = spacetime.db( p=p, DS="bigmemory.filenames" )
      rootdir = file.path( p$project.root, "interpolated" )
      dir.create( rootdir, showWarnings=FALSE, recursive =TRUE)
      fn.S =  file.path( rootdir, paste( "spacetime", "statistics", p$spatial.domain, "rdata", sep=".") )

      if ( DS=="statistics" ) {
        stats = NULL
        if (file.exists( fn.S) ) load( fn.S )
        return( stats )
      }

      if ( DS=="statistics.bigmemory.initialize" ) {
        # statistics storage matrix ( aggregation window, coords ) .. no inputs required
        spacetime.db( p=p, DS="bigmemory", grp="statistics.coordinates"  ) #Sloc
        spacetime.db( p=p, DS="bigmemory", grp="statistics.results"  ) # S
        return( "complete" )
      }

      if ( DS=="statistics.bigmemory.size" ) {
        # load bigmemory data objects pointers
        p = spacetime.db( p=p, DS="bigmemory.filenames" )
        S = bigmemory::attach.big.matrix(p$descriptorfile.S , path=p$tmp.datadir )
        return( nrow(S) )
      }

      if ( DS=="statistics.bigmemory.status" ) {
        # find locations for statistic computation and trim area based on availability of data
        # stats:
        p = spacetime.db( p=p, DS="bigmemory.filenames" )
        S = bigmemory::attach.big.matrix(p$descriptorfile.S , path=p$tmp.datadir )

        bnds = try( spacetime.db( p, DS="boundary" ) )

        if (!is.null(bnds)) {
          if( !("try-error" %in% class(bnds) ) ) {
            # problematic and/or no data (e.g., land, etc.) and skipped
            to.ignore =  which( bnds$inside.polygon == 0 ) # outside boundary
            i = which( is.nan( S[,1] ) & bnds$inside.polygon != 0 )
            # not yet completed
            j = which( is.na( S[,1] )  & bnds$inside.polygon != 0 )
            # completed
            k = which( is.finite (S[,1])  & bnds$inside.polygon != 0 ) # not yet done
        } } else {
            to.ignore = NA
            i = which( is.nan( S[,1] )  )
            # not yet completed
            j = which( is.na( S[,1] )   )
            # completed
            k = which( is.finite (S[,1])  ) # not yet done
        }

        return( list(problematic=i, incomplete=j, completed=k, n.total=nrow(S[]),
                     n.incomplete=length(j), n.problematic=length(i), n.complete=length(k), to.ignore=to.ignore ) )
      }

      if ( DS =="statistics.redo" ) {
        #\\ spacetime.db( "statsitics.redo") .. statistics are stored at a different resolution than the final grid
        #\\   this fast interpolates the solutions to the final grid
        p = spacetime.db( p=p, DS="bigmemory.filenames" )
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
          RES = spacetime.interpolate.kernel.density( x=p$sbbox$plons, y=p$sbbox$plats, z=z,
            locsout=locsout,  nxout=length(p$plons), nyout=length( p$plats),
            theta=p$dist.mwin, xwidth=p$dist.mwin*10, ywidth=p$dist.mwin*10 ) # 10 SD of the normal kernel
          # 10 SD of the normal kernel
          if ( !is.null( RES )) stats[,iv] = exp( RES$z ) # return to correct scale

          method = FALSE
          if (method=="inla.fast") { # fast, but not fast enough for prime time yet
            # interpolation using inla is also an option
            # but will require a little more tweaking as it was a bit slow
            range0 = median( ss$range, na.rm=TRUE )
            oo = which( is.finite( ss[,vn] ) )
            if ( length(oo) < 30 ) next()
            RES = spacetime.interpolate.inla.singlepass (
              Y=ss[oo,vn], locs=ss[oo, p$variables$LOCS], plocs=locsout, method="fast", link=datalink[iv] )
            if ( !is.null( RES )) stats[,iv] = RES$xmean
            rm (RES); gc()
          }
        }

        save( stats,  file=fn.S, compress=TRUE )
        return( fn.S)

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
    }

  }

