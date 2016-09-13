
  spacetime.db = function( DS, p, B=NULL, grp=NULL ) {
    #// usage: low level function to convert data into file-based data obects to permit parallel
    #// data access and manipulation and deletes/updates
    #// B is the xyz  or xytz data to work upon
    #// bigmemory operations are faster but do not support arrays .. left her in case read/write speed becomes an issue

    if (DS %in% "filenames" ) {
      # input data stored as a hdf5 file to permit operations with min memory usage
      # split into separate components to minimize filelocking conflicts
      p$tmp.datadir = file.path( p$project.root, "tmp" )
      if( !file.exists(p$tmp.datadir)) dir.create( p$tmp.datadir, recursive=TRUE, showWarnings=FALSE )
      p$ptr =list()
      p$ptr$Y =    file.path( p$tmp.datadir, "input.Y.hdf5")
      p$ptr$Ycov = file.path( p$tmp.datadir, "input.Ycov.hdf5" )
      p$ptr$Yloc = file.path( p$tmp.datadir, "input.Yloc.hdf5")
      p$ptr$P =    file.path( p$tmp.datadir, "predictions.hdf5")
      p$ptr$Pcov = file.path( p$tmp.datadir, "predictions_cov.hdf5")
      p$ptr$Ploc = file.path( p$tmp.datadir, "predictions_loc.hdf5")
      p$ptr$S =    file.path( p$tmp.datadir, "statistics.hdf5")
      p$ptr$Sloc = file.path( p$tmp.datadir, "statistics_loc.hdf5")
      return(p)
    }

    # --------------------------
    if (DS %in% "filelist" ) {
      p = spacetime.db( p=p, DS="filenames" )
      fl = c( 
           p$ptr$P,
           p$ptr$Ploc,
           p$ptr$Pcov,
           p$ptr$S,
           p$ptr$Sloc,
           p$ptr$Y,
           p$ptr$Ycov,
           p$ptr$Yloc
      )
      return( fl )
    }

    # --------------------------
    if (DS %in% "cleanup" ) {
      todelete = spacetime.db(p=p, DS="filelist")
      for (fn in todelete ) if (file.exists(fn)) file.remove(fn)
      return( todelete )
    }


    # ------------------
    if (DS == "data.initialize" ) {

      p = spacetime.db( p=p, DS="filenames" ) 

      if ( "data.frame" %in% class(B) ) {
        # B = B # nothing to do
      } else if ( "SpatialGridDataFrame" %in% class(B) ) {
        B = slot(B, "data") 
      }

      # dependent variable
      if ( file.exists( p$ptr$Y ) ) file.remove( p$ptr$Y )
      Y_ = B[, p$variables$Y ] 
      if ( exists("spacetime.link", p) ) Y_ = p$spacetime.link ( Y_ ) 
      Y = h5file(name =p$ptr$Y, mode = "a")
      Y["Y", compression=0L, chunksize=p$hdf5.chunksize] = as.vector(Y_)
      h5close(Y)
      rm(Y_)

      # independent variables/ covariate
      if (exists("COV", p$variables)) {
        if ( file.exists( p$ptr$Ycov ) ) file.remove( p$ptr$Ycov )
        Ycov = h5file(name =p$ptr$Ycov, mode = "a")
        Ycov["Ycov", compression=0L, chunksize=p$hdf5.chunksize] = as.matrix( B[ , p$variables$COV ] )
        h5close(Ycov)
      }

     # data coordinates
      if ( file.exists( p$ptr$Yloc ) ) file.remove( p$ptr$Yloc )
      Yloc = h5file(name =p$ptr$Yloc, mode = "a")
      Yloc["Yloc", compression=0L, chunksize=p$hdf5.chunksize] = as.matrix( B[ , p$variables$LOCS ])
      h5close(Yloc)
  
      return( "complete" )
    }


    #---------------------
    if (DS == "predictions.initialize" ) {

      p = spacetime.db( p=p, DS="filenames" ) 
      
      if ( "data.frame" %in% class(B) ) {
        # B = B # nothing to do
      } else if ( "SpatialGridDataFrame" %in% class(B) ) {
        B = slot(B, "data")  
      }
      
      # prediction coordinates
      if ( file.exists( p$ptr$Ploc ) ) file.remove( p$ptr$Ploc )
      Ploc = h5file(name =p$ptr$Ploc, mode = "a")
      Ploc["Ploc", compression=0L, chunksize=p$hdf5.chunksize] = as.matrix( B[, p$variables$LOCS ] )
      h5close(Ploc)

      # prediction covariates i.e., independent variables/ covariates
      if (exists("COV", p$variables)) {
        if ( file.exists( p$ptr$Pcov ) ) file.remove( p$ptr$Pcov )
        Pcov = h5file(name =p$ptr$Pcov, mode = "a")
        Pcov["Pcov", compression=0L, chunksize=p$hdf5.chunksize] = as.matrix( B[ , p$variables$COV ] )
        h5close(Pcov)
      }

      # predictions and associated stats
      if ( file.exists( p$ptr$P ) ) file.remove( p$ptr$P )
      P = h5file(name =p$ptr$P, mode = "a")
      P["P", compression=0L, chunksize=p$hdf5.chunksize] = matrix( NA_real_, nrow=p$nPreds, ncol=3 ) # pred, count, sd
      h5close(P)
      
      return( "complete" )
    }



    # -----------------

    if (DS %in% c( "boundary.redo", "boundary" ) )  {
      p = spacetime.db( p=p, DS="filenames" )
      fn =  file.path(p$tmp.datadir, "boundary.rdata" )
      if (DS=="boundary") {
        boundary = NULL
        if( file.exists(fn)) load( fn)
        return( boundary )
      }

      # data:
      Y = h5file( p$ptr$Y)["Y"]
      hasdata = 1:length(Y[])
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA
      h5close(Y)

      # covariates (independent vars)
      if ( exists( "COV", p$variables) ) {
        Ycov = h5file( p$ptr$Ycov)["Ycov"]
        if ( length( p$variables$COV ) == 1 ) {
          bad = which( !is.finite( Ycov) )
        } else {
          bad = which( !is.finite( rowSums(Ycov)) )
        }
        if (length(bad)> 0 ) hasdata[bad] = NA
        h5close(Ycov)
      }

      ii = na.omit(hasdata)
      ndata = length(ii)

      Yloc = h5file( p$ptr$Yloc)["Yloc"]
      # locs_noise = Yloc[ii,] + runif( ndata*2, min=-p$pres*p$spacetime.noise, max=p$pres*p$spacetime.noise )
      maxdist = max( diff( range( Yloc[ii,1] )), diff( range( Yloc[ii,2] )) )

      convex = -0.04
      if (exists( "mesh.boundary.convex", p) ) convex=p$mesh.boundary.convex
      resolution = 125
      if (exists( "mesh.boundary.resolution", p) ) resolution=p$mesh.boundary.resolution
      boundary=list( polygon = inla.nonconvex.hull(  Yloc[ii,], convex=convex, resolution=resolution ) )
      Sloc = h5file( p$ptr$Sloc)["Sloc"] # statistical output locations
      boundary$inside.polygon = point.in.polygon( Sloc[,1], Sloc[,2],
          boundary$polygon$loc[,1], boundary$polygon$loc[,2], mode.checked=TRUE)
      h5close(Sloc)
      save( boundary, file=fn, compress=TRUE )
      plot( Yloc[ii,], pch="." ) # data locations
      lines( boundary$polygon$loc , col="green" )
      h5close(Yloc)      
      return( fn )
    }


    # -----------------
    
    if (DS %in% c( "statistics.initialize", "statistics.size" , "statistics.status", "statistics.box" ) ) {
      p = spacetime.db( p=p, DS="filenames" )  
     
      if ( DS=="statistics.initialize" ) {
        # statistics storage matrix ( aggregation window, coords ) .. no inputs required
        # statistics coordinates
        Sloc_ = expand.grid( p$sbbox$plons, p$sbbox$plats )
        Sloc_ = as.matrix( Sloc_)
        if ( file.exists( p$ptr$Sloc ) ) file.remove( p$ptr$Sloc )
        Sloc = h5file(name =p$ptr$Sloc, mode = "a")
        Sloc["Sloc", compression=0L, chunksize=p$hdf5.chunksize] <- Sloc_
        h5close(Sloc)

        S_ = matrix( NA_real_, nrow=nrow(Sloc_), ncol=length( p$statsvars ) )
        if ( file.exists( p$ptr$S ) ) file.remove( p$ptr$S )
        S = h5file(name =p$ptr$S, mode = "a")
        S["S", compression=0L, chunksize=p$hdf5.chunksize] <- S_
        h5close(S)
        return( "complete" )
      }
      

      if ( DS=="statistics.size" ) {
        S = h5file( p$ptr$S)["S"]
        nS =  nrow(S) 
        h5close(S) 
        return(nS)
      }
      
      if (DS == "statistics.box")  {
        sbbox = list( plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$dist.mwin ),
                      plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$dist.mwin ) )
        return(sbbox)
      }


      if ( DS=="statistics.status" ) {
        # find locations for statistic computation and trim area based on availability of data
        # stats:
        S = h5file( p$ptr$S)["S"]
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

        out = list(problematic=i, incomplete=j, completed=k, n.total=nrow(S) ,
                     n.incomplete=length(j), n.problematic=length(i), 
                     n.complete=length(k), to.ignore=to.ignore )
        out$prop_incomp=out$n.incomplete / ( out$n.problematic + out$n.incomplete + out$n.complete)
        h5close(S)
        cat( paste("Proportion incomplete:", round(out$prop_incomp,5), "\n" )) 
        return( out )
      }
    }
  }

