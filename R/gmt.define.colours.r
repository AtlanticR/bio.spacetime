
  gmt.define.colours = function(params, variable="") {
    # "-I" inverts the colours
    # "-Z" continuous spectrum
      blue2red = c( "temp", "sal", "mr", "mrT", "smr", "smrT",
                     "mrPvalue", "mrPvalueT", "ca1", "ca2" )
      red2blue = NA
      if (variable %in% red2blue) {
        params$colourflags = "-Z"
      } else {
        params$colourflags = "-Z -I" # the I inverts the colour sequence
      }
      if (! is.null(params$gmtcol) ) {
        if (params$gmtcol == "red2green") params$colourflags = "-Z"
      }
    return (params)
  }


