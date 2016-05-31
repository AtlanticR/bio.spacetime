
  gstat.model.variogram = function(vgm.e, vp) {
    vgm.m=NULL
    all= NULL
    if (vp$skip) return(NULL)
    if (vp$fit) {
      vmodels = c( "Sph", "Exp", "Cir", "Pen" )
    }

    # initial start -- based upon empirical observations of stable variograms
    vgm.e$id = "var1"
    vp$psill = max(vgm.e$gamma) * 0.2
    vp$nugget = max(vgm.e$gamma) * 0.3
    vp$range = max(vgm.e$dist) *0.25

    vp0 = vp

    # find estimates of nugget and psill
    all = find.best.variogram.model( vmodels, vp, vgm.e )
    all = all[ which(all$range < (max(vgm.e$dist) ) ) ,] # likely solutions
    all = all[ which(all$range > (min(vgm.e$dist) ) ) ,] # likely solutions
    #      all = all[ which(all$psill < (max(vgm.e$gamma) ) ) ,] # likely solutions
    
    use.default = F 
    if (nrow(all) > 0 ) {
      good = all[ which.min(all$sse), ]
      vgm.m = vgm( psill=good$psill, model=good$model, range=good$range, nugget=good$nugget, kappa=good$kappa )
      vgm.m = fit.variogram(vgm.e, model=vgm.m)
      if ( vgm.m$range[2] > (max(vgm.e$dist) ) )  use.default=T
      if ( vgm.m$range[2] < (min(vgm.e$dist) ) )  use.default=T
    }
    
    if (is.null( vgm.m) | use.default ) { # still null .. use a simple default Spherical with conservative settings
      vgm.m = vgm( psill=vp0$psill, model="Sph", range=vp0$range, nugget=vp0$nugget )
      # vgm.m = fit.variogram(vgm.e, model=vgm.m)
    }
    return (vgm.m)
  }


