
spatial_warp = function( Z0, L0, L1, p0, p1, method="fast" ) {
  #\\ regrid/reproject from p0 to p1 ; rgdal calls this "warping" ;)

  if (method=="raster") {
    # here Z0 has to be rasterizable or a list of rasters or a brick
    for (vn in names(Z0)) {
       Z[[vn]] = raster::projectRaster(
          from = raster::rasterize( Z0, spatial_parameters_to_raster(p0), field=vn, fun=mean),
          to   = bio.spacetime::spatial_parameters_to_raster( p1) )
    }
    return (Z)
  }

  # ----------

  if (method=="fast") {
      # extract coords, convert to new coords, interpolate where required an
      M = matrix( NA, nrow=p0$nplons, ncol=p0$nplats) # matrix respresentation of the data in new coord system
      L0$plon = grid.internal( L0$plon, p0$plons ) # ensure correct resolution
      L0$plat = grid.internal( L0$plat, p0$plats )
      L2M = round( cbind( ( L0$plon-p0$plons[1])/p0$pres + 1, (L0$plat-p0$plats[1])/p0$pres + 1)) # row, col indices in matrix form of the new coordinate system
      L1 = planar2lonlat( L1, proj.type=p1$internal.projection )
      L1 = lonlat2planar( L1, proj.type=p0$internal.projection )  # convert lon, lat to old projection
      M[L2M] = Z0
      Z = fields::interp.surface( list( x=p0$plons, y=p0$plats, z=M ), loc=L1[,c("plon","plat")] )
      ii = which( !is.finite( Z ) )
      if ( length( ii) > 0 ) {
        # try again ..
        Z[ii] = fields::interp.surface( list( x=p0$plons, y=p0$plats, z=M ), loc=L1[ ii, c("plon","plat")] )
      }
      ii = which( !is.finite( Z ) )
      if ( length( ii) > 0 ) {
        if ( !exists( "wght", p1 ) ) {
          p1$wght = fields::setup.image.smooth(
            nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres,
            theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )
        }
        Zii =  fields::image.smooth( M, dx=p1$pres, dy=p1$pres, wght=p1$wght )$z
        Z[ii] = Zii[ii]
      }
      return( Z)
  }



if (0) {
  
        p1 = spatial_parameters( type=gr ) #target projection
           
        if ( p0$spatial.domain != p1$spatial.domain ) {

          Z = expand.grid( plon=p1$plons, plat=p1$plats, KEEP.OUT.ATTRS=FALSE )
          Zi = as.matrix( round( cbind( 
            ( Z$plon-p1$plons[1])/p1$pres + 1, (Z$plat-p1$plats[1])/p1$pres + 1
          ) ) ) 
     
          Z = planar2lonlat( Z, proj.type=p1$internal.crs )
          Z$plon_1 = Z$plon # store original coords
          Z$plat_1 = Z$plat
          Z = lonlat2planar( Z, proj.type=p0$internal.crs )
          p1_wgts = fields::setup.image.smooth( 
            nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres,
            theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )
           
          for (vn in varnames) {
            M = matrix(NA, nrow=p0$nplons, ncol=p0$nplats )
            M[Z0i] = Z0[[vn]]
            Znew = fields::interp.surface( list(x=p0$plons, y=p0$plats, z=M), loc=Z[, c("plon", "plat")] ) #linear interpolation
            Z[[vn]] = c(Znew)
            ii = which( !is.finite( Z[[vn]] ) )
            if ( length( ii) > 0 ) {
              MM = matrix(NA, nrow=p1$nplons, ncol=p1$nplats )
              MM[Zi] = Z[[vn]]
              Znew = fields::image.smooth( MM, dx=p1$pres, dy=p1$pres, wght=p1_wgts )
              Zii = fields::interp.surface( list(x=p1$plons, y=p1$plats, z=Znew$z), loc=Z[, c("plon_1", "plat_1")] ) #linear interpolation
              Z[[vn]][ii] = Zii[ii]
            }
          }
          Z$plon = Z$plon_1
          Z$plat = Z$plat_1
        
        } else {
          Z = Z0
        }

        Z$plon_1 = Z$plat_1 = Z$lon = Z$lat = NULL

}

}


