
  spacetime.db = function( DS, p, B=NULL, grp=NULL ) {
    #// usage: low level function to convert data into file-based data obects to permit parallel
    #// data access and manipulation and deletes/updates
    #// B is the xyz  or xytz data to work upon
    #// bigmemory operations are faster but do not support arrays .. left her in case read/write speed becomes an issue

    if (DS %in% "filenames" ) {
      # input data stored as a ff file to permit operations with min memory usage
      # split into separate components to reduce filelocking conflicts
      p$tmp.datadir = file.path( p$project.root, "tmp" )
      if( !file.exists(p$tmp.datadir)) dir.create( p$tmp.datadir, recursive=TRUE, showWarnings=FALSE )
      p$ptr =list()
      p$ptr$Y =    file.path( p$tmp.datadir, "input.Y.ff")
      p$ptr$Ycov = file.path( p$tmp.datadir, "input.Ycov.ff" )
      p$ptr$Yloc = file.path( p$tmp.datadir, "input.Yloc.ff")
      p$ptr$P =    file.path( p$tmp.datadir, "predictions.ff")
      p$ptr$Psd =  file.path( p$tmp.datadir, "predictions_sd.ff")
      p$ptr$Pn =   file.path( p$tmp.datadir, "predictions_n.ff")
      p$ptr$Pcov = file.path( p$tmp.datadir, "predictions_cov.ff")
      p$ptr$Ploc = file.path( p$tmp.datadir, "predictions_loc.ff")
      p$ptr$S =    file.path( p$tmp.datadir, "statistics.ff")
      p$ptr$Sloc = file.path( p$tmp.datadir, "statistics_loc.ff")
      return(p)
    }


    # --------------------------
    if (DS %in% "cleanup" ) {
      p = spacetime.db( p=p, DS="filenames" )
      for (fn in p$ptr ) if (file.exists(fn)) file.remove(fn)
      return( "done" )
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
      Y_ = B[, p$variables$Y ] 
      if ( exists("spacetime.link", p) ) Y_ = p$spacetime.link ( Y_ ) 
      Y = ff( Y_, filename=p$ptr$Y, overwrite=TRUE )
      close(Y) # ensure write/sync
      p$ff$Y = Y  # store pointer

      # independent variables/ covariate
      if (exists("COV", p$variables)) {
        Ycov_ = as.matrix( B[ , p$variables$COV ] )
        Ycov = ff::ff( Ycov_, dim=dim(Ycov_), filename=p$ptr$Ycov, overwrite=TRUE )
        close(Ycov)
        p$ff$Ycov = Ycov # store pointer
      }

     # data coordinates
      Yloc_ = as.matrix( B[ , p$variables$LOCS ])
      Yloc = ff::ff( Yloc_, dim=dim(Yloc_), filename=p$ptr$Yloc, overwrite=TRUE )
      close(Yloc)
      p$ff$Yloc = Yloc # store pointer
  
      return( p ) #return pointers to data
    }


    #---------------------
    if (DS == "predictions.initialize.xy" ) {

      p = spacetime.db( p=p, DS="filenames" ) 
      
      if ( "data.frame" %in% class(B) ) {
        # B = B # nothing to do
      } else if ( "SpatialGridDataFrame" %in% class(B) ) {
        B = slot(B, "data")  
      }
      
      # prediction coordinates
      Ploc_ = as.matrix( B[, p$variables$LOCS ] )
      Ploc = ff::ff( Ploc_, dim=dim(Ploc_), filename=p$ptr$Ploc, overwrite=TRUE )
      close(Ploc)
      p$ff$Ploc = Ploc  # store pointer

      # prediction covariates i.e., independent variables/ covariates
      if (exists("COV", p$variables)) {
        Pcov_ = as.matrix( B[ , p$variables$COV ] )
        Pcov = ff::ff( Pcov_, dim=dim(Pcov_), filename=p$ptr$Pcov, overwrite=TRUE )
        close(Pcov)
        p$ff$Pcov = Pcov  # store pointer
      }

      # predictions and associated stats
      P_ = matrix( NA, nrow=p$nPreds, ncol=3 ) # pred, count, sd
      P = ff::ff( P_, dim=dim(P), filename=p$ptr$P, overwrite=TRUE )
      close(P)
      p$ff$P = P  # store pointer
      
      return( p )
    }



    #---------------------
    if (DS == "predictions.initialize.xyt" ) {
      # prediction covariates i.e., independent variables/ covariates

      p = spacetime.db( p=p, DS="filenames" ) 

      if (exists("COV", p$variables)) {
        p$ff$Pcov = list()
        for ( i in p$variables$COV ){
          Pcov_ = as.matrix( B[[i]] )  # must pass as lists due to size constraints of ff
          fni = paste(p$ptr$Pcov, i, sep="_")
          p$ff$Pcov[i] = ff::ff( Pcov_, dim=dim(Pcov_), filename=fni, overwrite=TRUE )
          close(p$ff$Pcov[i])
        }
      }
      
      # predictions and associated stats
      # predictions and associated stats
      P_ = array( NA, dim=c(p$nplons, p$nplats, p$ny, p$nw) ) # 3=pred, count, sd
      P = ff::ff( P_, dim=dim(P_), filename=p$ptr$P, overwrite=TRUE )
      close(P)      
      p$ff$P = P  # store pointer

      Pn = ff::ff( P_, dim=dim(P_), filename=p$ptr$Pn, overwrite=TRUE )
      close(Pn)      
      p$ff$Pn = Pn  # store pointer
      
      Psd = ff::ff( P_, dim=dim(P_), filename=p$ptr$Psd, overwrite=TRUE )
      close(Psd)      
      p$ff$Psd = Psd  # store pointer
         
      return( p )
    }


    #---------------------
    if (DS == "predictions.initialize.xyts" ) {

      p = spacetime.db( p=p, DS="filenames" ) 
      
      if (exists("COV", p$variables)) {
        p$ff$Pcov = list()
        for ( i in p$variables$COV ){
          Pcov_ = as.matrix( B[[i]] ) # B must be a list that is already formatted properly
          fni = paste(p$ptr$Pcov, i, sep="_")
          p$ff$Pcov[i] = ff::ff( Pcov_, dim=dim(Pcov_), filename=fni, overwrite=TRUE )
          close(p$ff$Pcov[i])
        }
      }
      
      # predictions and associated stats
      P_ = array( NA, dim=c(p$nplons, p$nplats, p$ny, p$nw) ) # 3=pred, count, sd
      P = ff::ff( P_, dim=dim(P_), filename=p$ptr$P, overwrite=TRUE )
      close(P)      
      p$ff$P = P  # store pointer

      Pn = ff::ff( P_, dim=dim(P_), filename=p$ptr$Pn, overwrite=TRUE )
      close(Pn)      
      p$ff$Pn = Pn  # store pointer
      
      Psd = ff::ff( P_, dim=dim(P_), filename=p$ptr$Psd, overwrite=TRUE )
      close(Psd)      
      p$ff$Psd = Psd  # store pointer
            
      return( p )
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
      Y = p$ff$Y
      hasdata = 1:length(Y)
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA
      close(Y)

      # covariates (independent vars)
      if ( exists( "COV", p$variables) ) {
        Ycov = p$ff$Ycov 
        if ( length( p$variables$COV ) == 1 ) {
          bad = which( !is.finite( Ycov[]) )
        } else {
          bad = which( !is.finite( rowSums(Ycov[])) )
        }
        if (length(bad)> 0 ) hasdata[bad] = NA
        close(Ycov)
      }

      ii = na.omit(hasdata)
      ndata = length(ii)

      Yloc = p$ff$Yloc 
      # locs_noise = Yloc[ii,] + runif( ndata*2, min=-p$pres*p$spacetime.noise, max=p$pres*p$spacetime.noise )
      maxdist = max( diff( range( Yloc[ii,1] )), diff( range( Yloc[ii,2] )) )

      convex = -0.04
      if (exists( "mesh.boundary.convex", p) ) convex=p$mesh.boundary.convex
      resolution = 125
      if (exists( "mesh.boundary.resolution", p) ) resolution=p$mesh.boundary.resolution
      boundary=list( polygon = inla.nonconvex.hull(  Yloc[ii,], convex=convex, resolution=resolution ) )
      Sloc = p$ff$Sloc  # statistical output locations
      boundary$inside.polygon = point.in.polygon( Sloc[,1], Sloc[,2],
          boundary$polygon$loc[,1], boundary$polygon$loc[,2], mode.checked=TRUE)
      
      save( boundary, file=fn, compress=TRUE )
      plot( Yloc[ii,], pch="." ) # data locations
      lines( boundary$polygon$loc , col="green" )
      
      close(Sloc)
      close(Yloc)      
      return( fn )
    }


    # -----------------
    
    if (DS %in% c( "statistics.initialize", "statistics.size" , "statistics.status", "statistics.box" ) ) {
      p = spacetime.db( p=p, DS="filenames" )  
     
      if ( DS=="statistics.initialize" ) {
        # statistics storage matrix ( aggregation window, coords ) .. no inputs required
        # statistics coordinates
        Sloc_ = as.matrix( expand.grid( p$sbbox$plons, p$sbbox$plats ))
        Sloc = ff::ff( Sloc_, dim=dim(Sloc_), filename=p$ptr$Sloc_, overwrite=TRUE )
        close(Sloc)
        p$ff$Sloc = Sloc

        S_ = matrix( NA, nrow=nrow(Sloc_), ncol=length( p$statsvars ) )
        S = ff::ff( S_, dim=dim(S_), filename=p$ptr$S, overwrite=TRUE )
        close(S)
        p$ff$S = S

        return( p )
      }
      

      if ( DS=="statistics.size" ) {
        nS = nrow(p$ff$S) 
        close(S) 
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
        S = p$ff$S
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
        close(S)
        cat( paste("Proportion incomplete:", round(out$prop_incomp,5), "\n" )) 
        return( out )
      }
    }
  }

