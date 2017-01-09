
spatial_warp = function( Z0, L0, L1, p0, p1, method="fast", L0i=NULL, L1i=NULL ) {
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

    L0_mat = matrix(NA, nrow=p0$nplons, ncol=p0$nplats )
    if (is.null(L0i)) L0i = array_map( "xy->2", L0[, c("plon", "plat")],
      corner=c(p0$plons[1], p0$plats[1]), res=c(p0$pres, p0$pres) )
    L0_mat[L0i] = Z0
    L0_grid = list(x=p0$plons, y=p0$plats, z=L0_mat)

    L1_interp = fields::interp.surface( L0_grid, loc=L1[, c("plon", "plat")] ) #linear interpolation
    ii = which( !is.finite( L1_interp ) )
    if ( length( ii) > 0 ) {
      L1_mat = matrix(NA, nrow=p1$nplons, ncol=p1$nplats )
      if (is.null(L1i)) L1i = array_map( "xy->2", L1[, c("plon", "plat")],
        corner=c(p1$plons[1], p1$plats[1]), res=c(p1$pres, p1$pres) )
      L1_mat[L1i] = L1_interp
      if (!exists("wght", p1)) p1$wght = fields::setup.image.smooth( nrow=p1$nplons,
        ncol=p1$nplats, dx=p1$pres, dy=p1$pres, theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )
      L1_sm = fields::image.smooth( L1_mat, dx=p1$pres, dy=p1$pres, wght=p1$wght )
      L1_sm_interp = fields::interp.surface( list(x=p1$plons, y=p1$plats, z=L1_sm$z), loc=L1[, c("plon_1", "plat_1")] ) #linear interpolation from smoothed surface
      L1_interp[ii] = L1_sm_interp[ii]
    }
    return( L1_interp)

  }

}


