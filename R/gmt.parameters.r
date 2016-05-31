
  gmt.parameters = function(p=NULL) {

    # print( "Parameterisations for GMT ver 5.X.X " )
    
    if( is.null(p)) p=list()

    # default gmt-related parameters for Maritimes region of Canada
    p$gmtcol       = "seis"  # "red2green" is an alternate

    # interpolation setttings 
    p$maskres      = "-S25k" # to define area of influence to  accept after interpolation
    p$interpres    = "-n40k" #  to define area of influence for interpolation
    p$tension      = "-T0.4"  # 0.35+ for steep; 0.25 for smooth .. tension interpolation parameter
    p$bathy.tension = "-T0.75"  # tension parameter for GMT splines-in-tension (1=harmonic surface ... max,min not exceeded)
    p$bathy.maskres  = "-S40k"        # resolution
    p$bathy.contour= "-C100 -S4 -W0.25p" # contour lines every 100m with lines of various thinknesses

    p$block = T  # do a block ?
    p$blocktype = "mean"  # do a block mean ?

    # postscript generation options
    p$polygon.options = "-W0.8p"
    p$coast.options = "-Df -G220 -W0.5p"  # -Df = full resolution and colour/thickness
    p$scale.location = "-D4.75i/0.75i/2.6i/0.16ih" # alternate: "-D4.5i/0.8i/2.5i/.25ih"
    p$fnt = "18p"
    p$psresolution = "300"
    p$delta  = 1
    p$overlay = ""
    
    p$annot.lon0 = -63.7
    p$annot.lat0 = 47.25
    p$annot.fontsize = 32 # points
    p$annot.fontno = 0 # try  pstext  -L for possibles .. 0 = helvetica
    p$annot.fontcol = 0 # try  pstext  -L for possibles .. 0 = helvetica
    p$annot.angle = 0
    p$annot.justify = "LT"  # justification
    p$annot.text = ""
    
    p$delete.postscript = T
    p$redo.basemap = F

    basemap.location = project.datadirectory("bathymetry", "maps" )
    p$basemap = file.path(basemap.location, "basemap.default.ps" )  # default

    if ( exists( "spatial.domain", p) ) {
      
      if ( p$spatial.domain == "snowcrab" ) {
        # 2min grids
        p$res = "-I2m"  
        p$res.isobaths = "-I30s"  # slightly higher resolution for bathymetry to make it look nice
        p$isobaths = c( 50, 150, 250, 350, 450 ) ## in GMT
        p$bathy.zrange = "-Sa510/NaN -Sb1/NaN"
        p$basemap = file.path(basemap.location, "basemap.snowcrab.ps")
        p$annot = "-B2neSW"
        p$overlay = c( "cfanorth", "cfasouth", "cfa4x") # polygons of fishing areas 
        p$gmtproj = "-JL-61.5/45/46.5/43.5/6.5i"  # Lambert.conformal.conic.crab 
        p$region = "-R-66.4/-57.2/42.2/47.4"
      }
      
      if ( p$spatial.domain %in% c("SSE", "4vwx", "4vw", "4v", "4x") ) {
        # 2min grids
        p$res = "-I2m"  
        p$res.isobaths = "-I30s"
        p$isobaths = c( 50, 150, 250, 350, 450 ) ## in GMT
        p$bathy.zrange = "-Sa510/NaN -Sb1/NaN"
        p$basemap = file.path(basemap.location, "basemap.SSE.ps")
        p$annot = "-B2neSW"
        p$gmtproj = "-JL-61/45/46/44/6.5i"  # Lambert.conformal.conic 
        p$region = "-R-68/-56.5/41.5/47.5"
      }
      
      if ( p$spatial.domain %in% c("ecnasap", "ecnasap2") ) {
        # 2min grids
        p$res = "-I2m"  
        p$res.isobaths = "-I30s"
        p$isobaths = c(  200, 400, 600, 800 ) ## in GMT
        p$bathy.zrange = "-Sa1010/NaN -Sb1/NaN"
        p$basemap = file.path(basemap.location, "basemap.ecnasap.ps")
        p$annot = "-B5neSW"
        if( p$spatial.domain == "ecnasap" ) {
          p$gmtproj = "-JL-55/50/40/60/6.5i"
          p$region = "-R-72/-40/36.5/67.5"
        } 
        if( p$spatial.domain == "ecnasap2" ) {
          p$gmtproj =  "-JL-58/50/45/55/6.5i" # Lambert.conformal.conic
          p$region = "-R-73/-43/39/61"
        }
      }
 
      if ( p$spatial.domain == "canada.east" ) {
        # 15 arc sec ... high res for 
        p$res = "-I10m"
        p$res.bathy = "-I15s"  # for bathymetry ... high res preferred
        p$res.isobaths = "-I15s"  # isobaths to match bathy
        p$isobaths = c(  200, 400, 600, 800 ) ## in GMT
        p$bathy.zrange = "-Sa1010/NaN -Sb1/NaN"
        p$basemap = file.path(basemap.location, "basemap.canada.east.ps")
        p$annot = "-B10neSW"
        p$gmtproj = "-JQ-62/6.5i"  # Canada.east.PlateCarre-62
        p$region = "-R-72/-52/40/50"
      }
   
      if ( p$spatial.domain == "canada.east.highres" ) {
        # 15 arc sec ... high res for 
        p$res = "-I10m"
        p$res.bathy = "-I10s"  # for bathymetry ... high res preferred
        p$res.isobaths = "-I10s"  # isobaths to match bathy
        p$isobaths = c(  200, 400, 600, 800 ) ## in GMT
        p$bathy.zrange = "-Sa1010/NaN -Sb1/NaN"
        p$basemap = file.path(basemap.location, "basemap.canada.east.highres.ps")
        p$annot = "-B10neSW"
        p$gmtproj = "-JQ-62/6.5i"  # Canada.east.PlateCarre-62
        p$region = "-R-72/-52/40/50"
      }

    } else {
      stop( " 'spatial.domain' not found ... this needs to be defined for GMT to work ")
    }
    
    return(p)
  }
 
