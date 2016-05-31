   find.best.variogram.model = function (vmodels, vp, vgm.e, fix.nugget=F) {
      r = NULL
      for (model in vmodels) {
        vm0 = NULL
        if (substring(model,1,3) == "Mat") {
          vp$kappa = as.numeric(substring(model, 4, nchar(model)))
          vp$var.mod = "Mat"
          vm0 = vgm( psill=vp$psill, model=vp$var.mod, range=vp$range, nugget=vp$nugget, kappa=vp$kappa )
        } else {
          vp$kappa = NA
          vp$var.mod = model
          vm0 = vgm( psill=vp$psill, model=vp$var.mod, range=vp$range, nugget=vp$nugget )
        }
        vmf = NULL
        if (! fix.nugget) {
          vmf = try ( fit.variogram(object=vgm.e, model=vm0), silent=T )
        } else  {
          vmf = try ( fit.variogram(object=vgm.e, model=vm0, fit.sill=F), silent=T )
        }

# fit.method: default fitting method, used by gstat. The default method uses
# weights $N_h/h^2$ with $N_h$ the number of point pairs and
# $h$ the distance. This criterion is not supported by theory,
# but by practice.  For other values of 'fit.method', see table
# 4.2 in the gstat manual.
# fit.method=5 is REML .. slower but more accurate

        if (!( "try-error" %in% class(vmf) ) )   {
            r = rbind(r, data.frame(cbind(
                    sse    = attr(vmf,"SSErr"),
                    model  = vp$var.mod,
                    range  = vmf$range[2],
                    psill  = vmf$psill[2],
                    nugget = vmf$psill[1],
                    kappa  = vp$kappa )))
       }
     }

      if (!is.null(r)) {
        r = factor2number(r, c("sse", "range", "psill", "nugget", "kappa"))
        r = factor2character(r, "model")
      }

      return (r)
    }


