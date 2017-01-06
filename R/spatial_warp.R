
spatial_warp = function( Z0, L0, L1, p0, p1, method="fast" ) {
  #\\ regrid/reproject from p0 to p1 ; rgdal calls this "warping" ;)
message("deprecated")

  if (method=="raster") {
    # here Z0 has to be rasterizable or a list of rasters or a brick
    for (vn in names(Z0)) {
       Z_vn = raster::projectRaster(
          from = raster::rasterize( Z0, spatial_parameters_to_raster(p0), field=vn, fun=mean),
          to   = bio.spacetime::spatial_parameters_to_raster( p1) )
    }
    return (Z)
  }

  # ----------

  if (method=="fast") {
      # extract coords, convert to new coords, interpolate where required an
      M = matrix( NA, nrow=p0$nplons, ncol=p0$nplats) # matrix respresentation of the data in new coord system
      L0i = as.matrix( round( cbind(
          ( L0$plon-p0$plons[1])/p0$pres + 1, (L0$plat-p0$plats[1])/p0$pres + 1 ) ) )

      L1i = as.matrix( round( cbind(
        ( L1$plon-p1$plons[1])/p1$pres + 1, (L1$plat-p1$plats[1])/p1$pres + 1
      ) ) )

      L1 = planar2lonlat( L1, proj.type=p1$internal.crs )
      L1$plon_1 = L1$plon # store original coords
      L1$plat_1 = L1$plat
      L1 = lonlat2planar( L1, proj.type=p0$internal.crs )
      p1_wgts = fields::setup.image.smooth(
        nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres,
        theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )

      M[Z0i] = Z0_vn

      Znew = fields::interp.surface( list(x=p0$plons, y=p0$plats, z=M), loc=L1[, c("plon", "plat")] ) #linear interpolation
      L1_vn = c(Znew)
      ii = which( !is.finite( L1_vn ) )
      if ( length( ii) > 0 ) {
        MM = matrix(NA, nrow=p1$nplons, ncol=p1$nplats )
        MM[L1i] = L1_vn
        Znew = fields::image.smooth( MM, dx=p1$pres, dy=p1$pres, wght=p1_wgts )
        Zii = fields::interp.surface( list(x=p1$plons, y=p1$plats, z=Znew$z), loc=L1[, c("plon_1", "plat_1")] ) #linear interpolation
        L1_vn[ii] = Zii[ii]
      }

      return( L1 )
  }

}


