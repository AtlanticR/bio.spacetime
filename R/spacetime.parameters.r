
  spacetime.parameters = function( stp=list() ) {
    
    # the following parameters are for inside and outside ... do not make them exact multiples as this seems to make things hang ..
    if ( !exists("inla.mesh.max.edge", stp))  stp$inla.mesh.max.edge = c(  0.025,   0.04 )    # proportion of 2*stp$dist.max or equivalent: c(inside,outside) -- must be positive valued
    if ( !exists("inla.mesh.offset", stp))  stp$inla.mesh.offset   = c( - 0.025,  - 0.05 )   # how much to extend inside and outside of boundary: proportion of dist.max .. neg val = proportion
    if ( !exists("inla.mesh.cutoff", stp)) stp$inla.mesh.cutoff   = c( - 0.05,   - 0.5 )    ## min distance allowed between points: proportion of dist.max ; neg val = proportion

    if ( !exists("inla.mesh.hull.radius", stp)) stp$inla.mesh.hull.radius = c( -0.04, - 0.08 ) ## radius of boundary finding algorythm ; neg val = proportion

    if ( !exists("inla.mesh.hull.resolution", stp)) stp$inla.mesh.hull.resolution = 125  ## resolution for discretization to find boundary

    if ( !exists("spacetime.noise", stp)) stp$spacetime.noise = 0.001  # add a little noise to coordinates to prevent a race condition

    if ( !exists("inla.alpha", stp)) stp$inla.alpha = 2 # bessel function curviness
    if ( !exists("inla.nsamples", stp)) stp$inla.nsamples = 5000 # posterior similations 
    if ( !exists("predict.in.one.go", stp)) stp$predict.in.one.go = FALSE # use false, one go is very very slow and a resource expensive method
    if ( !exists("predict.quantiles", stp)) stp$predict.quantiles = c(0.025, 0.975 )  # posterior predictions robustified by trimming extreme values 
    
    if ( !exists("debug.file", stp)) stp$debug.file = file.path( ecomod.workdirectory, "inla.debug.out" )

   return(stp)
  
  }


