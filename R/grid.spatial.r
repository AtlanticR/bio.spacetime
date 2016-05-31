
    grid.spatial = function( p, R=NULL, grid.method=block.mean ) {
      
			fn = file.path( project.datadirectory("bathymetry"), p$resolution.ofdata+extent ) 
     
      if ( is.null( R ) ) {
        load( fn )
        return( R )
      }
      
      R = bathymetry.db( p, DS="z.lonlat.rawdata" ) 
      R$lon = grid.internal( R$lon, p$lons )
      R$lat = grid.internal( R$lat, p$lats )
      R = R[is.finite(rowSums(R)) ,]

      gc()
      Z = block.spatial ( xyz=R, function.block=grid.mean ) 
      Z = xyz2grid( Z, p$lons, p$lats)
      save( Z, file=fn, compress=T ) 
      
      return( "Completed" )
    }


