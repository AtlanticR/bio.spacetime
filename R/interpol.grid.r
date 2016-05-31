
  interpol.grid = function(xyz, locs=NULL, params, method="tps", getdata=F, outfile="surface.tmp") {

    out = NULL
    out = with (params, {
  
    if (method=="kernel.density") {

    }

    if (method=="tps") {
      outfile = make.random.string(outfile)
      Surf = make.random.string(".gmt.surface.interpolated")
      Surf.xy = make.random.string(".gmt.surface.interpolated.xy")
      xyz.data = make.random.string("xyz.tmp")
      xyz.blocked = make.random.string("tmp.blocked")

      write.table( xyz, file=xyz.data, quote=F, col.names=F, row.names=F)
      

# ------------
# depth sign is correct to this point but changes to negative somewhere below ... why? ------
# ------------


      cmd( "blockmean", xyz.data, region, res, ">",  xyz.blocked )
      cmd( "surface", xyz.blocked, region, res, tension, paste("-G", Surf, sep=""))

#     Md = ".gmt.data.mask"
#     cmd( "grdmask", xyz.data, region, res,"-NNaN/1/1", maskres, paste("-G", Md,sep=""))

      Ml = make.random.string(".gmt.data.landmask")
      cmd( "grdlandmask", region, res,"-N1/NaN/NaN/NaN/NaN -Di", paste("-G", Ml, sep=""))

      cmd( "grdmath", Ml, Surf, "MUL =", outfile )
      cmd( "grd2xyz", outfile, "-s >", Surf.xy )

      out = read.table( file=Surf.xy, sep="\t", header=F, strip.white=T )
      names(out) = c("lon", "lat", "z")

      remove.files( c( xyz.data, outfile, Surf, Surf.xy, xyz.blocked ) )
    }

    if (method=="nneighbor") {
      # this method uses weights and search radius: use "K" for great circle distances
      outfile = make.random.string(outfile)
      tmp.xyz = make.random.string("tmp.xyz")
      cmd ( "nearneighbor", xyz, "-bi -W", Rsearch, region, res, tension, paste("-G", outfile, sep="") )
      cmd( "grd2xyz", outfile, region, ">", tmp.xyz)
      out = read.table(tmp.xyz, header=F)
      names(out) = c("lon", "lat", "z")
      remove.files ( c( tmp.xyz, outfile ))
    }

    if (method=="freeform") cmd(xyz)

    if (method == "interp") {
      out = interp( x=xyz[,1], y=xyz[,2], z=xyz[,3], xo=locs[,1], yo=locs[,2], extrap=F, duplicate="mean")$z
      if (outfile != "surface.tmp") write.table( out, file=outfile, quote=F, col.names=F, row.names=F)
    }

    if (method == "krige.gstat") {
      require(gstat)

#        params$drange=150
#        params$nmax=20
#        params$psill=0.2
#        params$model = "Bes"
#        params$nugget = 0.05 )
#        params$pres = 1

      data = data.frame(xyz)
      names(data) = c("x", "y", "z")
      newdata = expand.grid(x=locs[,1], y=locs[,2])

      g = gstat(id="variable", formula=z~1, loc=~x+y, data=data, maxdist=drange, nmax=nmax )
      g.v = variogram(g, cutoff=drange, cressie=T)
      g.vgm = vgm(psill=psill, model=model, range=drange, nugget=nugget)
      f = fit.variogram(g.v, g.vgm, print.SSE=T)
      plot(g.v, model=f)
      g = gstat(g, id="variable", model=f)
      out = predict(object=g, newdata=newdata, block=c(pres,pres) )
      if (outfile != "surface.tmp") write.table( out, file=outfile, quote=F, col.names=F, row.names=F)
    }

    if (method == "inv.dist.gstat") {
      # inverse distance weighted interpolation .. power ==
      require(gstat)
      data = data.frame(xyz)
      names(data) = c("x", "y", "z")
      newdata = expand.grid(x=locs[,1], y=locs[,2])
      g = gstat( id="variable", formula=z~1, loc=~x+y, data=data, maxdist=drange, nmax=nmax, set=list(idp = power)  )
      out = predict(object=g, newdata=newdata)
      if (outfile != "surface.tmp") write.table( out, file=outfile, quote=F, col.names=F, row.names=F)
    }

    return (out)
    })
    return (out)
  }


