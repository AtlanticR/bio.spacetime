  
gmt.basemap = function( p ) {
  # create the initial layer to which are added GMT's postscript layers
  if(! exists("isobaths", p)) stop( "Please define p$isobaths in gmt.parameters")
     # override defaults:
 		 tmpdir  = tempdir()
     p$psxy.options = "-W0.001p -Sp0.005p -Glightslategray"  # -Sp == point
     ps.new  = "-K > "  # start a new PS and continue
     ps.append = "-O -K >> "  # add to prexisting PS and continue
     ps.finalize = "-O >> "
     bathy.mask = file.path(tmpdir, make.random.string("gmt.depth.mask") )
     gmt.annotfile = file.path(tmpdir, make.random.string("gmt.annot") )
     gmt.clip = file.path(tmpdir, make.random.string("gmt.clip") )
     gmt.isobathdata = file.path(tmpdir, make.random.string(".gmt.isobathdata") )
     gmt.depths = file.path(tmpdir, make.random.string(".gmt.depths") )
     gmt.depth.mask = file.path(tmpdir, make.random.string(".gmt.depth.mask") )
     gmt.data.landmask = file.path(tmpdir, make.random.string(".gmt.data.landmask") )
     
     cmd( "gmt psbasemap", p$region, p$gmtproj, p$annot, ps.new, p$basemap )
     for (d in p$isobaths ){
       isob = NULL
       isob = isobath.db( p=p, DS="gmt", depths=d )
       # isob$symbol =
       if (is.null(isob)) next()
       write.table( isob, file=gmt.isobathdata, quote=FALSE, col.names=FALSE, row.names=FALSE )
       cmd( "gmt psxy", gmt.isobathdata, p$region, p$gmtproj, p$psxy.options, ps.append, p$basemap )
     }
     cmd( "gmt pscoast", p$region, p$gmtproj, "-Df -W0.1p", ps.append, p$basemap )
     # finalize with  dummy text
     dummy.annotation = " " 
     dummy.annot.options = "-N -F+f12p,Helvetica,black+cTL"
     write( dummy.annotation, file=gmt.annotfile )
     cmd( "gmt pstext", gmt.annotfile, dummy.annot.options, p$region, p$gmtproj, ps.finalize, p$basemap )  
  }


