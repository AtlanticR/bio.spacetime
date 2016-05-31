
 gmt.colourscale = function(params, y, var="", NSTD=3) {

    signif.digits = 2
    YYY = y # copy

    trimoff10percent = quantile(y, probs=c(0.05, 0.95), na.rm=T)
    y = y[ which( y > trimoff10percent[[1]] &  y < trimoff10percent[[2]] ) ]
    if (length(y) < 30) y = YYY  # do not trim if there is little data
    rmean = mean( y, na.rm=T )
    rsd = NSTD * sd( y, na.rm=T )
    r0 = signif( rmean - rsd, signif.digits)
    r1 = signif( rmean + rsd, signif.digits)
    range = r1-r0
    crange=paste("-T", r0, "/", r1, "/", signif(range/10, signif.digits), sep="" )
    incscale = paste("-B", signif(range/4, 1), sep="")

    
    if ( !is.null(params$colourscale.quantiles) && params$colourscale.quantiles) {
      quants = quantile(y, probs=params$colourscale.quantiles.probs, na.rm=T)
      r0 = signif( quants[1], signif.digits)
      r1 = signif( quants[2], signif.digits)
      range = r1-r0
      crange=paste("-T", r0, "/", r1, "/", signif(range/10, signif.digits), sep="" )
      incscale = paste("-B", signif(range/4, 1), sep="")
    }
  
  # exceptions
 
    # snow crab specific 
      #if (var=="cpue") {
      #  crange="-T0/150/.1"
      #  incscale = "-B50"
      #}
      #if (var=="effort") {
      #  crange="-T0/450/.5"
      #  incscale = "-B150"
      #}
      #if (var=="landings") {
      #  crange="-T0/50/.1"
      #  incscale = "-B10"
      #}
      #if (var=="pd_anom") {
      #  crange="-T0.3/0.7/0.01"
      #  incscale = "-B.1"
      #}
      if (unlist(strsplit(var, ".", fixed=T))[1] == "sexratio") {
        crange="-T0/1/0.01"
        incscale = "-B.25"
      }
      if (unlist(strsplit(var, ".", fixed=T))[1] == "t") {
        crange="-T-0.5/8/0.01"
        incscale = "-B2"
      }
     
      # groundfish specific
      #if ( length(grep("totno", var)) == 1 |  length(grep("totwgt", var)) == 1   ){
      #    r.l = 1.5
      #    r.u = 4.5
      #    incscale = paste("-B", signif((r.u-r.l)/3 ,1), sep="")
      #    crange=paste("-T", r.l, "/", r.u, "/", (r.u-r.l)/10, sep="" )
      #  if ( length(grep("demersal", var)) == 1 ) {
      #    r.l = 1.5
      #    r.u = 4.5
      #    incscale = paste("-B", signif((r.u-r.l)/3 ,1), sep="")
      #    crange=paste("-T", r.l, "/", r.u, "/", (r.u-r.l)/10, sep="" )
      #  }
      #  if ( length(grep("pelagic", var)) == 1 ) {
      #    r.l = 0.5
      #    r.u = 3.5
      #    incscale = paste("-B", signif((r.u-r.l)/3 ,1), sep="")
      #    crange=paste("-T", r.l, "/", r.u, "/", (r.u-r.l)/10, sep="" )
      #  }
      #}
      
      if (length(grep("rmean", var)) == 1) {
        NSTD = 1.5  # override
        crange = paste("-T", -NSTD, "/", NSTD, "/", NSTD/4, sep="")
        incscale = "-B0.5"
      }
      if (length(grep("pmean", var)) == 1) {
        r.l = 0
        r.u = 1
        crange=paste("-T", r.l, "/", r.u, "/", 0.1, sep="" )
        incscale = "-B.25"
      }
      #if ( length(grep("mmean", var)) == 1)  {
      #  r.l = 1
      #  r.u = 4
      #  crange=paste("-T", r.l, "/", r.u, "/", (r.u-r.l)/10, sep="" )
      #  incscale = paste("-I -B", signif((r.u-r.l)/2 ,1), sep="")
      #  incscale = paste("-I -B", 0.5, sep="")
      #}
      #if ( length(grep("lmean", var)) == 1)  {
      #  r.l = 1
      #  r.u = 3
      #  crange=paste("-T", r.l, "/", r.u, "/", (r.u-r.l)/10, sep="" )
      #  incscale = paste("-I -B", signif((r.u-r.l)/2 ,1), sep="")
      #  incscale = paste("-I -B", 0.5, sep="")
      #}


 #     if ( length(grep("mr", var)) == 1 ) {
 #       r.l = 0
 #       r.u = 3.5
 #       crange=paste("-T", r.l, "/", r.u, "/", (r.u-r.l)/10, sep="" )
 #       incscale = paste("-I -B", signif((r.u-r.l)/5 ,2), sep="")
 #     }

 #     if ( length(grep("smr", var)) == 1 ) {
 #       r.l = 0
 #       r.u = 1.0
 #       crange=paste("-T", r.l, "/", r.u, "/", (r.u-r.l)/10, sep="" )
 #       incscale = paste("-I -B", signif((r.u-r.l)/5 ,2), sep="")
 #     }
# 
 #     if ( length(grep("ca", var)) == 1) {
 #       r.l = -2.5
 #       r.u = 2.5
 #       crange=paste("-T", r.l, "/", r.u, "/", 0.25, sep="" )
 #       incscale = "-B1"
 #     }

      if ( length(grep("rsquared", var)) == 1) {
        r.l = 0
        r.u = 1.0
        crange=paste("-T", r.l, "/", r.u, "/", 0.1, sep="" )
        incscale = "-B0.25"
      }
      #if ( length(grep("b1", var)) == 1) {
      #  r.l = -2
      #  r.u = 0
      #  crange=paste("-T", r.l, "/", r.u, "/", 0.1, sep="" )
      #  incscale = "-B0.4"
      #}


    params$crange=crange
    params$incscale=incscale

    return(params)

  }


