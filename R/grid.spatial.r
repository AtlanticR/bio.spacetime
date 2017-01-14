
    grid.spatial = function( p, R=NULL, grid.method=block.mean ) {

			fn = file.path( project.datadirectory("bio.bathymetry"), p$resolution.ofdata+extent )

      if ( is.null( R ) ) {
        load( fn )
        return( R )
      }


 lons = seq(p$lon0, p$lon1, by=p$dres)
 lats = seq(p$lat0, p$lat1, by=p$dres)



      R = bathymetry.db( p, DS="z.lonlat.rawdata" )
      R$lon = grid.internal( R$lon, lons )
      R$lat = grid.internal( R$lat, lats )
      R = R[is.finite(rowSums(R)) ,]

      gc()
      Z = block.spatial ( xyz=R, function.block=grid.mean )
      Z = xyz2grid( Z, lons, lats)
      save( Z, file=fn, compress=T )

      return( "Completed" )
    }


