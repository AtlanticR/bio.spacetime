
  xyz2grid = function (xyz, xx, yy, fill=NA, method="R", params=NULL)  {

    if (method == "R") {

      # this function expects evenly distributed x,y locations (data must be blocked already)
      # this is a direct copy of a function from library fields (?)
      nr = nrow(xyz )
      nc = ncol(xyz )
      x = xyz[, 1]
      y = xyz[, 2]
      if ( nc==2 ) {
        z = rep.int(1, length.out=nr ) #add a dummy value
      } else {
        z = xyz[, 3]
      }
      rm(xyz); gc()
      nx = length(xx)
      ny = length(yy)
      xmin = min(xx)
      xmax = max(xx)
      xrange = xmax - xmin
      ymin = min(yy)
      ymax = max(yy)
      yrange = ymax - ymin
      dx = mean( diff(xx), na.rm=T)
      dy = mean( diff(yy), na.rm=T)
      out = matrix(fill, nrow=nx, ncol=ny)
      row = round((x - xmin)/dx) + 1
      col = round((y - ymin)/dy) + 1
      row.col = cbind(row, col)
      out[row.col] = z
    }

    if (method=="gmt") {
		
			tmpdir  = tempdir()
 
# not used ... and not complete ... do not use until checked
      require(ncvar)
      with(params, {
      # primarily intended for lon/lat but can be used for planar coords too .. check GMT options

      tmp.grid = file.path(tmpdir, make.random.string("tmp.grid") )
      tmp.xyz = file.path(tmpdir, make.random.string("tmp.xyz") )
      tmpfile =  file.path(tmpdir, make.random.string("tmp.file") )

#      write.table(xyz, file=tmpfile, quote=F, col.names=F, row.names=F )
#      cmd( "xyz2grd", tmp.xyz, paste("-G", tmp.grid, sep=""), region, res)
#      cmd( "grd2xyz", tmpfile, "-ZBa > xyz2.tmp")

#      dims =  var.get.ncv(path=tmpfile, name="dimension")
      #z = var.get.ncv(path=tmpfile, name="z")
#      z = read.table("xyz2.tmp", header=F)
#      out = matrix( data=z[,1], nrow=dims$data[1], ncol=dims$data[2] )

#      remove.files ( c( tmpfile, tmp.xyz ) )

      }) # end with params
    }

    if (method=="gmt1") {

      # primarily intended for lon/lat but can be used for planar coords too .. check GMT options
		
			tmpdir  = tempdir()
 
      require(ncvar)
      out = with(params, {
        xyz.tmp = file.path(tmpdir, make.random.string("xyz.tmp"))
        tmp.grid = file.path(tmpdir, make.random.string("tmp.grid"))
        tmp2.grid = file.path(tmpdir, make.random.string("tmp2.grid"))
        write.table(xyz, file=xyz.tmp, quote=F, col.names=F, row.names=F )
        #cmd( "xyz2grd", xyz.tmp, paste("-G", tmp.grid, sep=""), region, res)
#        cmd( "blockmode", xyz.tmp, xyz2.tmp, region, res )
        cmd( "surface", xyz.tmp, paste("-G", tmp.grid, sep=""), region, res, "-T0.35" )
        cmd( "grdsample", tmp.grid, paste("-G", tmp2.grid, sep=""),  region, res)
        out = grid.netcdf(tmp2.grid)
        remove.files ( c( xyz.tmp, tmp.grid, tmp2.grid ) )
        return(out)
      }) # end with
    }

    return (out)
    }


