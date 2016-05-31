

spatial.parameters = function( p=NULL, type=NULL ) {

  if (is.null(p)) p=list()

  if ( ! exists("spatial.domain", p) ) p$spatial.domain = type
  if ( ! is.null(type)) p$spatial.domain = type  # type has priority over p$spatial.domain
   
  if ( p$spatial.domain %in% c("SSE", "snowcrab") ) {
    # source raw data for bathymetry:
		p$bathymetry.xyz = file.path( project.datadirectory("bathymetry"), "data", "bathymetry.canada.east.xyz" )  # ascii
		p$bathymetry.bin = file.path( project.datadirectory("bathymetry"), "data", "bathymetry.canada.east.bin" )  # GMT binary
		# resolution and region
		p$internal.projection = "utm20"
    p$internal.crs =  "+proj=utm +ellps=WGS84 +zone=20 +units=km"
    p$dres = 1/60/4  # this is the 15 second grid from CHS  .. default use highest resolution
    p$pres = 1   # 1 km resolution!
    p$lon0=-68
    p$lon1=-56
    p$lat0=41
    p$lat1=48
    p$lons = seq(p$lon0, p$lon1, by=p$dres)
    p$lats = seq(p$lat0, p$lat1, by=p$dres)
    p$nlons = length(p$lons)
    p$nlats = length(p$lats)
    p$corners = data.frame(lon=c(p$lon0,p$lon1), lat=c(p$lat0,p$lat1))
    p$corners = lonlat2planar( p$corners, proj.type=p$internal.projection )
    
    p$corners$plon = round( p$corners$plon, 0)  # this matches the p$pres value of 1 km resolution
    p$corners$plat = round( p$corners$plat, 0)  # this matches the p$pres value of 1 km resolution
    
    p$plons = seq(min(p$corners$plon), max(p$corners$plon), by=p$pres)
    p$plats = seq(min(p$corners$plat), max(p$corners$plat), by=p$pres)
    p$nplons = length(p$plons)
    p$nplats = length(p$plats)
  }
  
  if ( p$spatial.domain %in% c("SSE.mpa") ) {
    # source raw data for bathymetry:
		p$bathymetry.xyz = file.path( project.datadirectory("bathymetry"), "data", "bathymetry.canada.east.xyz" )  # ascii
		p$bathymetry.bin = file.path( project.datadirectory("bathymetry"), "data", "bathymetry.canada.east.bin" )  # GMT binary
		# resolution and region
		p$internal.projection = "utm20"
    p$internal.crs =  "+proj=utm +ellps=WGS84 +zone=20 +units=km"
    p$dres = 1/60/4  # this is the 15 second grid from CHS  .. default use highest resolution
    p$pres = 1  # 1 km resolution
    p$lon0=-71
    p$lon1=-48
    p$lat0=37
    p$lat1=48
    p$lons = seq(p$lon0, p$lon1, by=p$dres)
    p$lats = seq(p$lat0, p$lat1, by=p$dres)
    p$nlons = length(p$lons)
    p$nlats = length(p$lats)
    p$corners = data.frame(lon=c(p$lon0,p$lon1), lat=c(p$lat0,p$lat1))
    p$corners = lonlat2planar( p$corners, proj.type=p$internal.projection )
    
    p$corners$plon = round( p$corners$plon, 0)  # this matches the p$pres value of 1 km resolution
    p$corners$plat = round( p$corners$plat, 0)  # this matches the p$pres value of 1 km resolution
    
    p$plons = seq(min(p$corners$plon), max(p$corners$plon), by=p$pres)
    p$plats = seq(min(p$corners$plat), max(p$corners$plat), by=p$pres)
    p$nplons = length(p$plons)
    p$nplats = length(p$plats)
    p$boundary = "scotia.fundy.with.buffer"
  }

  if ( p$spatial.domain=="canada.east") {
		# source raw data for bathymetry:
		p$bathymetry.xyz = file.path( project.datadirectory("bathymetry"), "data", "bathymetry.canada.east.xyz" )  # ascii
		p$bathymetry.bin = file.path( project.datadirectory("bathymetry"), "data", "bathymetry.canada.east.bin" )  # GMT binary
		# resolution and region
		p$internal.projection = "lambert.conic.canada.east"
    p$internal.crs = "+proj=lcc +ellps=WGS84  +lon_0=62W +lat_0=45N +lat_1=43N +lat_2=47N +units=km "
    p$dres = 1/60/4  # this is the 15 second grid from CHS  .. ~ 0.25 km
    p$pres = 1  # km
    p$lon0=-72
    p$lon1=-52
    p$lat0=40
    p$lat1=50
    p$lons = seq(p$lon0, p$lon1, by=p$dres)
    p$lats = seq(p$lat0, p$lat1, by=p$dres)
    p$nlons = length(p$lons)
    p$nlats = length(p$lats)
    p$corners = data.frame(lon=c(p$lon0,p$lon1), lat=c(p$lat0,p$lat1))
    p$corners = lonlat2planar( p$corners, proj.type=p$internal.projection )
      
    p$corners$plon = round( p$corners$plon, 0)  # this matches the p$pres value of 1 km resolution
    p$corners$plat = round( p$corners$plat, 0)  # this matches the p$pres value of 1 km resolution
    
    p$plons = seq(min(p$corners$plon), max(p$corners$plon), by=p$pres)
    p$plats = seq(min(p$corners$plat), max(p$corners$plat), by=p$pres)
    p$nplons = length(p$plons)
    p$nplats = length(p$plats)
  }

  if ( p$spatial.domain %in% c("canada.east.highres")) {
		# source raw data for bathymetry:
		p$bathymetry.xyz = file.path( project.datadirectory("bathymetry"), "data", "bathymetry.canada.east.xyz" )  # ascii
		p$bathymetry.bin = file.path( project.datadirectory("bathymetry"), "data", "bathymetry.canada.east.bin" )  # GMT binary
		# resolution and region
		p$internal.projection = "lambert.conic.canada.east"
    p$internal.crs = "+proj=lcc +ellps=WGS84  +lon_0=62W +lat_0=45N +lat_1=43N +lat_2=47N +units=km"
    p$dres = 1/60/4  # CHS is 15 arc second ~ 0.25 km
    p$pres = 0.5  # discretize to 0.5 km resolution
    p$lon0=-72
    p$lon1=-52
    p$lat0=40
    p$lat1=50
    p$lons = seq(p$lon0, p$lon1, by=p$dres)
    p$lats = seq(p$lat0, p$lat1, by=p$dres)
    p$nlons = length(p$lons)
    p$nlats = length(p$lats)
    p$corners = data.frame(lon=c(p$lon0,p$lon1), lat=c(p$lat0,p$lat1))
    p$corners = lonlat2planar( p$corners, proj.type=p$internal.projection ) 
    
    # this must be sufficient to capture the p$pres value of 0.25 km resolution
    p$corners$plon = round( p$corners$plon, 1) # 1 because the resolution is 0.5 km 
    p$corners$plat = round( p$corners$plat, 1)  
    
    p$plons = seq(min(p$corners$plon), max(p$corners$plon), by=p$pres)
    p$plats = seq(min(p$corners$plat), max(p$corners$plat), by=p$pres)
    p$nplons = length(p$plons)
    p$nplats = length(p$plats)
    
    
  }

  return(p)
}


