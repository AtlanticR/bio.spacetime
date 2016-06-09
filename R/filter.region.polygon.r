
  filter.region.polygon = function( x, region, planar=F, proj.type="utm20" ) {

    region = recode.areas( region )
    out = NULL
    for (reg in region) {
     poly = read.table( find.bio.gis(reg), header=F)
      names(poly) =c("lon", "lat")

      a = NULL

      if (planar) {
        poly.planar = lonlat2planar (poly, proj.type=proj.type)
        a = which(point.in.polygon(x$plon, x$plat, poly.planar$plon, poly.planar$plat) != 0 )
      }

      if (!planar) {
        a = which(point.in.polygon(x$lon, x$lat, poly$lon, poly$lat) != 0 )
      }

      out = c(out, a)
    }
    out = sort(unique(out))
    return(out)
  }


