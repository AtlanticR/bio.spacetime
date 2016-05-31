  lookup.projection.params = function(x) {
    # for examples, see: 
    # http://www.progonos.com/furuti/MapProj/Normal/TOC/cartTOC.html
    out  = switch( x,
      utm20             = "+proj=utm +ellps=WGS84 +zone=20 +units=km ", 
      utm20.substrate   = "+proj=utm +datum=NAD83 +zone=20 +units=km",
      lambert.conic.equidist    = "+proj=eqdc +ellps=WGS84 +lon_0=63W +lat_0=45N +lat_1=43N +lat_2=47N +units=km ",
      lambert.conic             = "+proj=lcc +ellps=WGS84  +lon_0=63W +lat_0=45N +lat_1=43N +lat_2=47N +units=km ",
      lambert.conic.4t          = "+proj=lcc +ellps=WGS84  +lon_0=63W +lat_0=47N +lat_1=45N +lat_2=49N +units=km ",
      lambert.conic.canada.east = "+proj=lcc +ellps=WGS84  +lon_0=62W +lat_0=45N +lat_1=43N +lat_2=47N +units=km ",
      tmercator         = "+proj=tmerc +ellps=WGS84 +x_0=500000 +y=-400000 +lon_0=90w +units=km "
    )
    return ( out )
  }
