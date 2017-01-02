
spatial_parameters = function( p=NULL, type=NULL ) {

    if (is.null(p)) p=list()

    if ( ! exists("spatial.domain", p) ) p$spatial.domain = type
    if ( ! is.null(type)) p$spatial.domain = type  # type has priority over p$spatial.domain

    if ( p$spatial.domain %in% c("SSE", "snowcrab") ) {
        p$internal.projection = "utm20"
        p$internal.crs =  "+proj=utm +ellps=WGS84 +zone=20 +units=km"
        p$dres = 1/60/4  # this is the 15 second grid from CHS  .. default use highest resolution
        p$pres = 1   # 1 km resolution!
        p$lon0=-68
        p$lon1=-56
        p$lat0=41
        p$lat1=48
        p$psignif = 1
    }

    if ( p$spatial.domain %in% c("SSE.mpa") ) {
        p$internal.projection = "utm20"
        p$internal.crs =  "+proj=utm +ellps=WGS84 +zone=20 +units=km"
        p$dres = 1/60/4  # this is the 15 second grid from CHS  .. default use highest resolution
        p$pres = 1  # 1 km resolution
        p$lon0=-71
        p$lon1=-48
        p$lat0=37
        p$lat1=48
        p$psignif = 1
        p$boundary = "scotia.fundy.with.buffer"
    }

    if ( p$spatial.domain=="canada.east") {
        p$internal.projection = "lambert.conic.canada.east"
        p$internal.crs = "+proj=lcc +ellps=WGS84  +lon_0=62W +lat_0=45N +lat_1=43N +lat_2=47N +units=km "
        p$dres = 1/60/4  # this is the 15 second grid from CHS  .. ~ 0.25 km
        p$pres = 1  # km
        p$lon0=-72
        p$lon1=-52
        p$lat0=40
        p$lat1=50
        p$psignif = 1
    }

    if ( p$spatial.domain %in% c("canada.east.highres")) {
        p$internal.projection = "lambert.conic.canada.east"
        p$internal.crs = "+proj=lcc +ellps=WGS84  +lon_0=62W +lat_0=45N +lat_1=43N +lat_2=47N +units=km"
        p$dres = 1/60/4  # CHS is 15 arc second ~ 0.25 km
        p$pres = 0.5  # discretize to 0.5 km resolution
        p$lon0=-72
        p$lon1=-52
        p$lat0=40
        p$lat1=50
        p$psignif = 1
    }

    if ( p$spatial.domain %in% c("canada.east.superhighres")) {
        p$internal.projection = "lambert.conic.canada.east"
        p$internal.crs = "+proj=lcc +ellps=WGS84  +lon_0=62W +lat_0=45N +lat_1=43N +lat_2=47N +units=km"
        p$dres = 1/60/4  # CHS is 15 arc second ~ 0.25 km
        p$pres = 0.2  # discretize to 0.2 km resolution
        p$lon0=-72
        p$lon1=-52
        p$lat0=40
        p$lat1=50
        p$psignif = 1
    }

    p$lons = seq(p$lon0, p$lon1, by=p$dres)
    p$lats = seq(p$lat0, p$lat1, by=p$dres)
    p$nlons = length(p$lons)
    p$nlats = length(p$lats)
    p$corners = data.frame(lon=c(p$lon0,p$lon1), lat=c(p$lat0,p$lat1))
    p$corners = lonlat2planar( p$corners, proj.type=p$internal.projection )
    p$corners$plon = signif( p$corners$plon, p$psignif)  # this matches the p$pres value of x km resolution
    p$corners$plat = signif( p$corners$plat, p$psignif)  # this matches the p$pres value of x km resolution
    p$plons = seq(min(p$corners$plon), max(p$corners$plon), by=p$pres)
    p$plats = seq(min(p$corners$plat), max(p$corners$plat), by=p$pres)
    p$nplons = length(p$plons)
    p$nplats = length(p$plats)

  return(p)
}


