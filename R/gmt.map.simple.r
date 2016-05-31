
 
  gmt.map.simple = function( gmt ) {

    dir.create ( dirname(gmt$out), recursive=T, showWarnings=F )
    
    gmt.prep (fnt="18p", psres=300 )
    cmd( "gmtset PAGE_ORIENTATION portrait")

    # tmp local GMT files
    cpt   = tempfile("gmt")  
    grid   = tempfile("gmt") 
    bin  = tempfile("gmt")  
    clip   = tempfile("gmt")
    block   = tempfile("gmt")
    land.mask   = tempfile("gmt")
    annot   = tempfile("gmt")
    depth.mask =  tempfile("gmt")
    mask =  tempfile("gmt")
    psout = gmt$out
    bathy.bin = tempfile("gmt")  
    bathy.block = tempfile("gmt")  
    bathy.mask = tempfile("gmt")  

 
    # define some local parameters
    append = "-O -K"

    # data check
    if ( !is.null(gmt$dat) ) {
      gmt$inp = tempfile("gmt")
      gmt$dat = gmt$dat[ is.finite(rowSums(gmt$dat)) ,]
      write.table( gmt$dat, file=gmt$inp, quote=F, col.names=F, row.names=F)
    } 

    # colour scale
    cmd( "makecpt", gmt$cpt, ">", cpt  )
   
    # discretise and interpolate
    cmd( "blockmedian", gmt$region, gmt$resolution, gmt$inp, ">", block )
    cmd( "surface", gmt$region, gmt$resolution, block, gmt$tension, paste("-G", grid, sep="") )

    # begin mapping
    cmd( "psbasemap", gmt$region, gmt$gmtproj, "-Bnesw -K >", psout )
    cmd( "psmask ",  gmt$region, gmt$gmtproj, gmt$resolution, block, append, "-S10K >>", psout ) # define which grids are reliable/not reliable

      if ( "linecontour" %in% gmt$outputs ) {
        cmd( "grdcontour", grid, gmt$gmtproj, paste( "-C", cpt, sep=""), append, "-A- -S0.2m -W+faint -L-350/-1 >>", psout )
      }
      if ( "colourcontour"  %in% gmt$outputs ) {
        cmd( "grdimage", grid, gmt$gmtproj, paste( "-C", cpt, sep=""), append, "-Q -nb/0.5 >>", psout )
      } 
      
    cmd( "psmask -C", append, ">>", psout )

    if ( "coastline" %in% gmt$outputs ) {
      cmd( "pscoast", gmt$region, gmt$gmtproj, append, "-Df -W0.1p >>", psout )
    }
    
    if ( "colourscale" %in% gmt$outputs ) {
      cmd( "psscale", paste("-C", cpt, sep=""), gmt$incscale, gmt$scale.location, append, ">>", psout ) # colourscale
    }
    
    if ( "annotation" %in% gmt$outputs ) {
      annotation = paste( gmt$annot.base, gmt$annot.text )
      write( annotation, file=annot )
      cmd( "pstext", annot, gmt$region, gmt$gmtproj, append, ">>", psout )  # variablename, year
    }
   
     # finalize with  dummy text
     dummy.annotation = " " 
     dummy.annot.options = "-N -F+f12p,Helvetica,black+cTL"
     write( dummy.annotation, file=annot )
     cmd( "gmt pstext", annot, dummy.annot.options, gmt$region, gmt$gmtproj, "-O >>", psout )  
     
    for (i in  c( cpt, grid, bin, clip, block, annot )) {
      rm(i)
    }
    
    print( gmt$out )

   }


