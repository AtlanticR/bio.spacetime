
  gmt.xyplot = function( params, data, year, conversions="ps2png" ) {

		tmpdir  = tempdir()
 
    with(params, {

    # tmp local GMT files
      cpt   = file.path(tmpdir, make.random.string(".gmt.cpt"))  # palette file
      grid  = file.path(tmpdir, make.random.string("gmt.grid"))
      bin   = file.path(tmpdir, make.random.string(".gmt.bin"))  # binary data
      clip  = file.path(tmpdir, make.random.string(".gmt.clip"))
      indat = file.path(tmpdir, make.random.string(".gmt.indat"))
      blocked = file.path(tmpdir, make.random.string(".gmt.blocked"))
      annotations = file.path(tmpdir, make.random.string(".gmt.annotations") )
      outfile = paste( file.path(tmpdir, make.random.string(".gmt.outfile")), ".ps", sep="") # final postscript file

      dir.create(basedir, recursive=T, showWarnings=F)
      write.table( data, file=indat, quote=F, col.names=F, row.names=F)

    # define some local parameters
      append = "-O -K"
      gmt.prep (fnt=fnt, psres=psresolution)

    # begin mapping
      cmd( "psbasemap", region, gmtproj, annot, "-K >", outfile )

    # overlay data
      cmd( "cat", basemap, ">>", outfile ) # bathymetry
      if ( length(overlay)>0 & !is.null(overlay)) {
      for (o in overlay) {
        cmd( "psxy", find.ecomod.gis(o), region, gmtproj, append, polygon.options, ">>", outfile ) # polygons
      }}
      cmd( "pscoast", region, gmtproj, coast.options, append, ">>", outfile ) # coastline
      cmd( "psxy", indat,  region, gmtproj, psxyoptions, append, " >>", outfile ) # data
      if (annot.text=="") annot.text = year
      if (as.character( annot.text) == "1000" ) annot.text = "Mean state"
      write( paste(annot.lon0, annot.lat0, annot.fontsize, annot.angle, annot.fontno, annot.justify, annot.text), file=annotations )
      cmd( "pstext", annotations, region, gmtproj, "-O >>", outfile )  # variablename, year
      
      remove.files ( c( cpt, grid, bin, clip, indat, blocked, annotations ) )
        
      if (!is.finite(conversions)) { #assume postscript not needed
        graphics.convert (outfile, outfile.basename, conversions)
      } else {  # assume just the postscript file is required
        file.copy(from=outfile, to=paste(outfile.basename,"ps",sep="."), overwrite=T)
      }

    })
      
    gmt.cleanup ()

    return ()
  }


