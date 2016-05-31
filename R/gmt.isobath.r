
  gmt.isobath = function( params, depth, redo.bin=F ) {
  # obtain co-ordinates of a given isobath
    params = gmt.parameters(params)  # get defaults ..
		tmpdir  = tempdir()
    a = with(params, {
      append = "-O -K"
      basemap = file.path(tmpdir, make.random.string("gmt.isobath.basemap")) 
        # not needed but retained to store default output from grdcontour
      gmt.bin = gsub(".xyz$", ".bin", inp)
      gmt.clip = file.path(tmpdir, make.random.string(".gmt.clip"))
      gmt.depths = file.path(tmpdir, make.random.string(".gmt.depths"))
      gmt.depth.mask = file.path(tmpdir, make.random.string(".gmt.depth.mask"))
      gmt.data.landmask = file.path(tmpdir, make.random.string(".gmt.data.landmask"))
      bathy.mask = file.path(tmpdir, make.random.string("gmt.depth.mask"))
      tmp.iso = file.path(tmpdir, make.random.string("tmp.iso"))
      isobaths = file.path(tmpdir, make.random.string("tmp.isobaths"))
      bathy.contour= paste( "-m -S4", paste("-D", isobaths, sep=""), paste("-C",depth, sep="") ) # override default
      bathy.zrange = paste("-Sa-1/NaN -Sb-", depth+1, "/NaN", sep="")
      # make the basemap with bathymetry 
      if (redo.bin) cmd( "gmtconvert -bo", inp, ">", gmt.bin )
      cmd( "blockmedian -bi3 -bo", gmt.bin, region, res, ">", gmt.clip )
      cmd( "surface -bi3", gmt.clip, region, res, bathy.tension, paste("-G", gmt.depths, sep="" ))
      cmd( "psmask -bi3", gmt.bin, region, res, bathy.maskres, gmtproj, append, ">", basemap )
      cmd( "grdclip", gmt.depths, bathy.zrange, paste("-G", gmt.depth.mask, sep="") )
      cmd( "grdlandmask", region, res,"-N1/NaN/NaN/NaN/NaN -Di", paste("-G", gmt.data.landmask, sep=""))
      cmd( "grdmath", gmt.data.landmask, gmt.depth.mask, "MUL =", bathy.mask)
      cmd( "grdcontour", bathy.mask, gmtproj, bathy.contour, append, ">>", basemap )
      cmd( "gawk '!/>/' ", isobaths, ">", tmp.iso )
      A = read.table(tmp.iso)
      remove.files ( c( basemap, gmt.bin, gmt.clip, gmt.depths, gmt.depth.mask, gmt.data.landmask, bathy.mask, tmp.iso ) )
      return(A)
    }) # end with

    return (a)
  }


