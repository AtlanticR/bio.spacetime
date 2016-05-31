 
  gmt.map.variables.core = function( id=NULL, U, params, variables, plottimes, basedir, conversions, delta, init.files, db ) {
   
    if (!is.null( init.files ) )  for (ii in init.files) source(ii) 
    varnames = colnames(U)
     
    # the first index is a list that is passed from the calling prog: in this case "ssplt" (if parallel)
    if (params$do.parallel) {
      id = as.numeric(id)
    } else {
      id = c(1:length(variables) )
    }
   
   
    for (i in id) {
      plotvar = variables[i]
      if (! (plotvar %in% varnames) ) next ()
      
      for (ti in plottimes) {
        u = recode.time(U, ti, delta )  # delta is for running averages
        u[,plotvar] = variable.recode( u[,plotvar], plotvar, "forward", db=db) # check if log transform is needed
        u = u[, c("yr", "lon", "lat", plotvar, "sa")]
        u = u[is.finite(u[,plotvar] *u[,"lon"]*u[,"lat"] ),]
        params$outdir = file.path(basedir, ti, plotvar)
        params = gmt.define.colours (params, variable=plotvar)
        params = gmt.colourscale(params, u[,plotvar], plotvar, NSTD=3 ) # NSTD is no of stdev
        dir.create ( params$outdir, recursive=T, showWarnings=F )
        for (j in sort(unique(u$yr)) ) {
          oo = which(u$yr==j & is.finite( u[,1] + u[,2] + u[,3]+u[,4] ) )
          if (length(oo) > 10 ) {
            toplot = u[ oo,  c("lon", "lat", plotvar, "sa") ]
            params$outfile.basename = file.path(params$outdir, paste(plotvar, j, sep="."))
            # print( params$outfile.basename )
            gmt.map( params, toplot, year=j, vname=plotvar, conversions=conversions )
          }
        }
    }}
    return ()
  }


