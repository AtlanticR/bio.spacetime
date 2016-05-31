
spacetime.regrid = function( Z0, L0, L1, p0, p1, method="fast" ) {
  #\\ regrid/reproject from p0 to p1 
  #\\ rgdal calls this "warping"
  
  if (method=="raster") {
    # here Z0 has to be rasterizable or a list of rasters or a brick
    require (raster)
    for (vn in names(Z0)) {
       Z[[vn]] = projectRaster( 
          from = rasterize( Z0, spatial.parameters.to.raster(p0), field=vn, fun=mean), 
          to   = spatial.parameters.to.raster( p1) )
    }
    return (Z)
  }

  if (method=="fast") {
      # extract coords, convert to new coords, interpolate where required an
      M = matrix( NA, nrow=p0$nplons, ncol=p0$nplats) # matrix respresentation of the data in new coord system
      L0$plon = grid.internal( L0$plon, p0$plons ) # ensure correct resolution
      L0$plat = grid.internal( L0$plat, p0$plats )
      L2M = cbind( ( L0$plon-p0$plons[1])/p0$pres + 1, (L0$plat-p0$plats[1])/p0$pres + 1) # row, col indices in matrix form of the new coordinate system
      L1 = planar2lonlat( L1, proj.type=p1$internal.projection )
      L1 = lonlat2planar( L1, proj.type=p0$internal.projection )  # convert lon, lat to old projection     
      M[L2M] = Z0
      Z = fields::interp.surface( list( x=p0$plons, y=p0$plats, z=M ), loc=L1[,c("plon","plat")] )
      ii = which( is.na( Z ) )
      if ( length( ii) > 0 ) {
        # try again ..
        Z[ii] = fields::interp.surface( list( x=p0$plons, y=p0$plats, z=M ), loc=L1[ ii, c("plon","plat")] )
      }
      ii = which( is.na( Z ) )
      if ( length( ii) > 0 ) {
        Zii =  fields::image.smooth( M, dx=p0$pres, dy=p0$pres, wght=p0$wght )$z  
        Z[ii] = Zii[ii]
      }
      return( Z)
  }
}


