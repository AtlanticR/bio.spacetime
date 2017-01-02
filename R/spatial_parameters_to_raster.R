
spatial_parameters_to_raster = function( params ) {
  #\\ Take a spatial parameter list wirh corners and resolution and CRS
  #\\ and convert to a raster template
  #\\ bio uses left edge as coordinates, raster uses center

  # if (edge.reference) {
  #   params$corners$plon = params$corners$plon + c(-0.5, 0.5)*params$pres
  #   params$corners$plat = params$corners$plat + c(-0.5, 0.5)*params$pres

  #   if (params$spatial.domain=="canada.east") {
  #     params$corners$plon = params$corners$plon + c(0, -0.5)*params$pres
  #     params$corners$plat = params$corners$plat + c(0, -0.5)*params$pres
  #   }
  #   if (params$spatial.domain=="canada.east.highres") {
  #     params$corners$plon = params$corners$plon + c(-1/2, -1/2)*params$pres
  #     params$corners$plat = params$corners$plat + c(-1/2, -1/2)*params$pres
  #   }
  #   if (params$spatial.domain=="canada.east.superhighres") {
  #     params$corners$plon = params$corners$plon + c(-1/2, -1/2)*params$pres
  #     params$corners$plat = params$corners$plat + c(-1/2, +1/2)*params$pres
  #   }
     
  # }
#  if  ( params$spatial.domain=="canada.east.superhighres" ) {
    params$corners$plon = params$corners$plon + c(-1/2, 0)*params$pres
    params$corners$plat = params$corners$plat + c(-1/2, +1/2)*params$pres # 
 # }

  ras = raster::raster(
    ncols=params$nplons,
    nrows=params$nplats,
    res=params$pres ,
    xmn= params$corners$plon[1], # rasters are center referenced
    xmx= params$corners$plon[2],
    ymn= params$corners$plat[1],
    ymx= params$corners$plat[2],
    # ext=extent ( rbind( params$corners$plon, params$corners$plat ) ),
    crs=params$internal.crs )


 # if ( (dim(ras)[1] != params$nplats) | (dim(ras)[2] != params$nplons) ) {
    print( dim(ras))
    print ( paste( params$nplats, params$nplons) )
    print( ras)
    print( params$corners )
 #   stop( "Parameters to raster is dimensionally not correct")
 # } 

  return(ras)

  if (0) {
    bioLibrary("bio.spacetime")
    require("rgdal")
    require("raster")


    params = spatial_parameters( type="canada.east.superhighres" )
    u = spatial_parameters_to_raster(params)

    params = spatial_parameters( type="canada.east.highres" )
    u = spatial_parameters_to_raster(params)

    params = spatial_parameters( type="canada.east" )
    u = spatial_parameters_to_raster(params)

    params = spatial_parameters( type="SSE" )
    u = spatial_parameters_to_raster(params)

    params = spatial_parameters( type="SSE.mpa" )
    u = spatial_parameters_to_raster(params)

    params = spatial_parameters( type="snowcrab" )
    u = spatial_parameters_to_raster(params)
  }

}


