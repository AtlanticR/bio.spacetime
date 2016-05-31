
spatial.parameters.to.raster = function( params, edge.reference=TRUE ) {
  #\\ Take a spatial parameter list wirh corners and resolution and CRS
  #\\ and convert to a raster template
  #\\ bio uses left edge as coordinates, raster uses center
  require( raster)
  if (edge.reference) {
    params$corners$plon = params$corners$plon + c(-0.5, 0.5)*params$pres
    params$corners$plat = params$corners$plat + c(-0.5, 0.5)*params$pres
    if (params$spatial.domain=="canada.east.highres") {
      params$corners$plon = params$corners$plon + c(0, -0.5)*params$pres
      params$corners$plat = params$corners$plat + c(0, -0.5)*params$pres
    }
  }

  ras = raster(
    ncols=params$nplons,
    nrows=params$nplats,
    res=params$pres ,
    xmn= params$corners$plon[1], # rasters are center referenced
    xmx= params$corners$plon[2],
    ymn= params$corners$plat[1],
    ymx= params$corners$plat[2],
    # ext=extent ( rbind( params$corners$plon, params$corners$plat ) ),
    crs=params$internal.crs )

    dd = dim(ras)
    if( dd[1] != params$nplats) stop( "Dim of plats in error")
    if( dd[2] != params$nplons) stop( "Dim of plons in error")

    ras_coord = coordinates(ras)
    uu = sort(unique(ras_coord[,1]))
    vv = sort(unique(ras_coord[,2]))

    u = unique( params$plons-uu )
    v = unique( params$plats-vv )

    up=FALSE # updated flag
    if ( length(u) > 1) stop( "problem 1")
    if ( length(u)==1 ) {
      if ( u!=0) {
        params$corners$plon = params$corners$plon + u
        up=TRUE
      }
    }

    if ( length(v) > 1) stop( "problem 2")
    if(length(v)==1){
      if ( v!=0) {
        params$corners$plat = params$corners$plat + v
        up=TRUE
    }}

    if(up) {
      print( "Raster coords offset issue ... trying a simple fix")
      print(paste(u, v))
      ras = raster(
        ncols=params$nplons,
        nrows=params$nplats,
        res=params$pres ,
        xmn= params$corners$plon[1], # rasters are center referenced
        xmx= params$corners$plon[2],
        ymn= params$corners$plat[1],
        ymx= params$corners$plat[2],
        # ext=extent ( rbind( params$corners$plon, params$corners$plat ) ),
        crs=params$internal.crs )

      # test again
      ras_coord = coordinates(ras)
      uu = sort(unique(ras_coord[,1]))
      vv = sort(unique(ras_coord[,2]))
      u = unique( params$plons-uu )
      v = unique( params$plats-vv )

      if( u!=0) {
        print(" This has been tested with known projections in bio, this is a new one? and raster library is not happy? You might need to add a new condition to spatial.parameters.toraster.. " )
        print(u)
        stop(" plons offsets are not unique/correct")
      }

      if( v!=0) {
        print(" This has been tested with known projections in bio, this is a new one? and raster library is not happy? You might need to add a new condition to spatial.parameters.toraster.." )
        print(v)
        stop(" plats offsets are not unique/correct? ")
      }
      print("...OK")
    }

  return(ras)

  if (0) {
    bioLibrary("bio.spacetime")
    require("rgdal")
    p = spatial.parameters( type="canada.east.highres" )
    u = spatial.parameters.to.raster(p)

    p = spatial.parameters( type="canada.east" )
    u = spatial.parameters.to.raster(p)

    p = spatial.parameters( type="SSE" )
    u = spatial.parameters.to.raster(p)

    p = spatial.parameters( type="SSE.mpa" )
    u = spatial.parameters.to.raster(p)

    p = spatial.parameters( type="snowcrab" )
    u = spatial.parameters.to.raster(p)
  }

}


