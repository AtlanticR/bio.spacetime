
  spacetime.db = function( DS, p, B=NULL, grp=NULL ) {
    #// usage: low level function to convert data into bigmemory obects to permit parallel
    #// data access and manipulation
    #// B is the xyz  or xytz data to work upon

    if (DS %in% "bigmemory.filenames" ) {
      # input data stored as a bigmatrix to permit operations with min memory usage
      # split into separate components to minimize filelocking conflicts
      p$tmp.datadir = file.path( p$project.root, "tmp" )
      if( !file.exists(p$tmp.datadir)) dir.create( p$tmp.datadir, recursive=TRUE, showWarnings=FALSE )

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

      p$backingfile.P_st = "predictions_st.bigmatrix.tmp"
      p$descriptorfile.P_st = "predictions_st.bigmatrix.desc"

      p$backingfile.Pcov_st = "predictions_cov_st.bigmatrix.tmp"
      p$descriptorfile.Pcov_st = "predictions_cov_st.bigmatrix.desc"

      return(p)
    }

    # --------------------------
    if (DS %in% "bigmemory.filelist" ) {
      p = spacetime.db( p=p, DS="bigmemory.filenames" )
      fl = file.path( p$tmp.datadir,
        c( p$backingfile.P, p$descriptorfile.P,
           p$backingfile.P_st, p$descriptorfile.P_st,
           p$backingfile.S, p$descriptorfile.S,
           p$backingfile.Sloc, p$descriptorfile.Sloc,
           p$backingfile.Ploc, p$descriptorfile.Ploc,
           p$backingfile.Pcov, p$descriptorfile.Pcov,
           p$backingfile.Pcov_st, p$descriptorfile.Pcov_st,
           p$backingfile.Y, p$descriptorfile.Y,
           p$backingfile.X, p$descriptorfile.X,
           p$backingfile.LOCS, p$descriptorfile.LOCS
      ))
      return( fl )
    }

    # --------------------------
    if (DS %in% "bigmemory.cleanup" ) {
      todelete = spacetime.db(p=p, DS="bigmemory.filelist")
      for (fn in todelete ) if (file.exists(fn)) file.remove(fn)
      return( todelete )
    }

    # ------------------
    if (DS == "bigmemory.dependent" ) {
      p = spacetime.db( p=p, DS="bigmemory.filenames" ) 
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
      return( fn.Y )
    }

    # ------------------
    if (DS == "bigmemory.covariates" ) {
      p = spacetime.db( p=p, DS="bigmemory.filenames" ) 
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
      return( fn.X )
    }

    # ------------------
    if (DS == "bigmemory.coordinates" ) {
      p = spacetime.db( p=p, DS="bigmemory.filenames" ) 
      # data coordinates
      fn.LOC = file.path(p$tmp.datadir, p$backingfile.LOC )
      if ( file.exists( fn.LOC) ) file.remove( fn.LOC)
      LOCS = filebacked.big.matrix( nrow=nrow(B), ncol=2, type="double", dimnames=NULL, separated=FALSE,
          backingpath=p$tmp.datadir, backingfile=p$backingfile.LOCS, descriptorfile=p$descriptorfile.LOCS )
      if ( "data.frame" %in% class(B) ) {
        LOCS[] = as.matrix( B[ , p$variables$LOCS ] )
      } else if ( "SpatialGridDataFrame" %in% class(B) ) {
        LOCS[] = as.matrix( coordinates(B) )
      }
      return( fn.LOC )
    }

    # ------------------

    if (DS == "bigmemory.prediction.coordinates" ) {
      p = spacetime.db( p=p, DS="bigmemory.filenames" ) 
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
      return( fn.Ploc )
    }

    # ------------------

    if (DS == "bigmemory.prediction.covariates" ) {
      p = spacetime.db( p=p, DS="bigmemory.filenames" ) 
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
      return( fn.Pcov )
    }

    # ------------------

    if (DS == "bigmemory.statistics.coordinates" ) {
      p = spacetime.db( p=p, DS="bigmemory.filenames" ) 
      # statistics coordinates
      fn.Sloc = file.path(p$tmp.datadir, p$backingfile.Sloc )
      if ( file.exists( fn.Sloc) ) file.remove( fn.Sloc )
      coords = expand.grid( p$sbbox$plons, p$sbbox$plats )
      Sloc = filebacked.big.matrix( nrow=nrow(coords), ncol=2, type="double", dimnames=NULL, separated=FALSE,
         backingpath=p$tmp.datadir, backingfile=p$backingfile.Sloc, descriptorfile=p$descriptorfile.Sloc )
      Sloc[] = as.matrix( coords )
      return( fn.Sloc )
    }

    # ------------------

    if (DS == "bigmemory.statistics.results" ) {
      p = spacetime.db( p=p, DS="bigmemory.filenames" ) 
      # statistics results output file .. initialize
      coords = expand.grid( p$sbbox$plons, p$sbbox$plats )
      fn.S = file.path(p$tmp.datadir, p$backingfile.S )
      if ( file.exists( fn.S) ) file.remove( fn.S)
      S = filebacked.big.matrix( nrow=nrow(coords), ncol= length( p$statsvars ), type="double", init=NA, dimnames=NULL, separated=FALSE,
        backingpath=p$tmp.datadir, backingfile=p$backingfile.S, descriptorfile=p$descriptorfile.S )
      return( fn.S )
    }
      
    # ----------------

    if (DS == "bigmemory.predictions"  ) {
      p = spacetime.db( p=p, DS="bigmemory.filenames" )
      rootdir = file.path( p$project.root, "interpolated" )
      dir.create( rootdir, showWarnings=FALSE, recursive =TRUE)
      # fn.P =  file.path( rootdir, paste( "spacetime", "predictions", p$spatial.domain, "rdata", sep=".") )
      # predictions storage matrix (discretized)
      fn.P = file.path(p$tmp.datadir, p$backingfile.P )
      if ( file.exists( fn.P) ) file.remove( fn.P)
      # contains c(count, pred.mean, pred.sd)
      P = filebacked.big.matrix( nrow=p$nplon * p$nplat, ncol=3, type="double", init=NA, dimnames=NULL, separated=FALSE,
        backingpath=p$tmp.datadir, backingfile=p$backingfile.P, descriptorfile=p$descriptorfile.P )
      return( fn.P )
    }

    # ----------------

    if (DS == "bigmemory.predictions_st"  ) {
      p = spacetime.db( p=p, DS="bigmemory.filenames" )
      rootdir = file.path( p$project.root, "interpolated" )
      dir.create( rootdir, showWarnings=FALSE, recursive =TRUE)
      # fn.P =  file.path( rootdir, paste( "spacetime", "predictions", p$spatial.domain, "rdata", sep=".") )
      # predictions storage matrix (discretized)
      fn.P = file.path(p$tmp.datadir, p$backingfile.P_st )
      if ( file.exists( fn.P_st) ) file.remove( fn.P_st)
      # contains c(count, pred.mean, pred.sd)
      P_st = filebacked.big.matrix( nrow=p$nplon * p$nplat * p$nyrs * p$ndyrs, ncol=3, type="double", init=NA, dimnames=NULL, separated=FALSE,
        backingpath=p$tmp.datadir, backingfile=p$backingfile.P_st, descriptorfile=p$descriptorfile.P_st )
      return( fn.P_st )
    }


    if (DS == "bigmemory.prediction.covariates_st" ) {
      p = spacetime.db( p=p, DS="bigmemory.filenames" ) 
      # prediction covariates i.e., independent variables/ covariates
      if ( exists( "X", p$variables) ) {
        fn.Pcov_st = file.path(p$tmp.datadir, p$backingfile.Pcov_st )
        if ( file.exists( fn.Pcov_st) ) file.remove( fn.Pcov_st)
        Pcov_st = filebacked.big.matrix( nrow=nrow(B), ncol=length( p$variables$X ), type="double", dimnames=NULL, separated=FALSE,
          backingpath=p$tmp.datadir, backingfile=p$backingfile.Pcov_st, descriptorfile=p$descriptorfile.Pcov_st )
        if ( "data.frame" %in% class(B) ) {
          Pcov_st[] = as.matrix( B[ , p$variables$X ] )
        } else if ( "SpatialGridDataFrame" %in% class(B) ) {
          Pcov_st[] = as.matrix( slot(B, "data")[, p$variables$X ]  )
        }
      }
      return( fn.Pcov )
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
    if (DS %in% c( "bigmemory.statistics", "bigmemory.statistics.size" , "bigmemory.statistics.status" ) ) {
      p = spacetime.db( p=p, DS="bigmemory.filenames" )
      rootdir = file.path( p$project.root, "interpolated" )
      dir.create( rootdir, showWarnings=FALSE, recursive =TRUE)
      fn.S =  file.path( rootdir, paste( "spacetime", "statistics", p$spatial.domain, "rdata", sep=".") )

      if ( DS=="bigmemory.statistics" ) {
        # statistics storage matrix ( aggregation window, coords ) .. no inputs required
        spacetime.db( p=p, DS="bigmemory", grp="statistics.coordinates"  ) #Sloc
        spacetime.db( p=p, DS="bigmemory", grp="statistics.results"  ) # S
        return( "complete" )
      }

      if ( DS=="bigmemory.statistics.size" ) {
        S = bigmemory::attach.big.matrix(p$descriptorfile.S , path=p$tmp.datadir )
        return( nrow(S) )
      }

      if ( DS=="bigmemory.statistics.status" ) {
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

    }

  }

