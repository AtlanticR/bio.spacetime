
  #\\ good references:
  #\\ http://www.stats.uwo.ca/faculty/aim/tsar/tsar.pdf
  #\\ http://wwwuser.gwdg.de/~cscherb1/content/Statistics%20Course%20files/A%20short%20introduction%20to%20time%20series%20analysis%20in%20R.pdf
  #\\ http://www.statoek.wiso.uni-goettingen.de/veranstaltungen/zeitreihen/sommer03/ts_r_intro.pdf
  #\\ http://www.uoguelph.ca/~ackerman/Primer-Ackerman-2010.pdf
  #\\ http://www.statsref.com/HTML/index.html?temporal_autocorrelation.html


  bioLibrary( "bio.spacetime" )

  #sunspot data
  o = timeseries.simulated.data( DS="sunspots.seasonal" )
  o = timeseries.simulated.data( DS="sunspots" )  # just annual

  o = timeseries.simulated.data( )  # random
  o$y = o$y0


# compare methods with sunspots data
  o = timeseries.simulated.data( DS="sunspots" )  # just annual
  sunspots = o$y
  time.diff = mean(diff(o$timestamp)) # unit of time
  plot(sunspots,type="b")

  x = sunspot.year
  x = sunspots
  o = timeseries.autocorrelation( x=x, freq=frequency(x) ) # default = spec.pgram
  o = timeseries.autocorrelation( x=x, method="fft", freq=frequency(x) )




