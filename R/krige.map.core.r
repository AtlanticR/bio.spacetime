krige.map.core = function( ip=NULL, p=NULL, init.files=NULL  ) {
   
  for (i in init.files) source( i )
  if ( is.null(id)) ip = c(1: p$nruns ) 

  for (iip in ip) {
        y = p$runs[iip,"y"]
        r = p$runs[iip,"r"]
        v = p$runs[iip,"v"]

    k.pname = paste(v, "pred", sep=".")
    k.vname = paste(v, "var",  sep=".")

#  these results have already been interpolated
#  .. reduce the degree of interpolation in the mapping of the results

    p$tension = "-T1"  # 0.35+ for steep; 0.25 for smooth
    p$maskres = "-S1k"  # effectively turns off interpolation
    p$interpres = "-n1k"

    S = get.PS.S.gridded (p, y, v)$S
    mest = variable.recode( S[,"kv"], v, direction="forward" )  # used for generating a colour scale that matches the whole time span
    rm (S); gc()
    
    bk0 = krigng.db( DS="UK.point.PS", p=list(v=v, y=y, r=r) )

    if (is.na(bk0) ) next ()
    bk = bk0$bk
    rm (bk0); gc()

    mdata = planar2lonlat( bk, proj.type=p$internal.projection ) # for mapping
#    mdata[,k.pname] = variable.recode( mdata[,k.pname], v, direction="backward")
#    mdata[,k.vname] = variable.recode( sqrt(mdata[,k.vname]), v, direction="backward")
    mdata[,k.vname] = sqrt(mdata[,k.vname])

    if (!file.exists(p$basemap) )  gmt.basemap(p)

# map means
    p$outdir = file.path( project.datadirectory("snowcrab"), "R", "predictions", "kriged.estimates")
    dir.create (path=p$outdir, recursive=T, showWarnings=F)
    p$outfile.basename = file.path(p$outdir, paste(v, y, sep="."))
    p = gmt.define.colours (p, v)
    p = gmt.colourscale(p, mest, v, NSTD=3 ) # NSTD is no of stdev .. uniform colour scale
    gmt.map( p, mdata[, c("lon", "lat", k.pname)], y, v , conversions="ps2png" )

# map variances
    p$outdir = file.path( project.datadirectory("snowcrab"), "R", "predictions", "kriged.variances")
    dir.create (path=p$outdir, recursive=T, showWarnings=F)
    p$outfile.basename = file.path(p$outdir, paste(v, y, sep="."))
    p = gmt.define.colours (p, v)
    p = gmt.colourscale(p, mdata[, k.vname], v, NSTD=3 ) # no need to force uniform scale as this is a diagnostic plot
    gmt.map( p, mdata[, c("lon", "lat", k.vname)], y, v, conversions="ps2png" )
   
  }
  return("Complete map")
}


