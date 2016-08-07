
require(sp)
require(microbenchmark)

data(meuse)
coordsK = data.frame( plon=meuse$y, plat=meuse$x ) 
dKK = as.matrix(dist( coordsK, diag=TRUE, upper=TRUE)) # distance matrix between knots

 # testing speed of MVN functions:
  sigma = c(1, 1)
  kappa = 1 
  phi =1
  rhoKs = exp(-phi * dKK)^kappa   ## spatial correlation
  covKs = sigma[2]*sigma[2] * rhoKs

  nKs = nrow(dKK)
  muKs = rnorm(nKs)
  muKs2 = rnorm(nKs)
  muKK = muKs %*% t(muKs2)

  microbenchmark( 
    LaplacesDemonCpp::rnorm(10e6, muKs, muKs2), 
    stats::rnorm(10e6, muKs, muKs2), 
    times=10 )

                                         expr      min       lq     mean   median       uq      max
 LaplacesDemonCpp::rnorm(1e+07, muKs, muKs2) 103.0105 103.5515 118.6936 105.3423 107.1295 176.2662
            stats::rnorm(1e+07, muKs, muKs2) 470.3615 476.6918 484.5429 478.1062 484.7872 536.9856
 
  
 microbenchmark( 
    LaplacesDemonCpp::dnorm_rcpp(muKs2,  muKs, log=TRUE), 
    stats::dnorm(muKs2,  muKs, log=TRUE), 
    times=10000 )

                                             expr    min      lq     mean  median     uq    max
 LaplacesDemonCpp::dnorm_rcpp(muKs2, muKs, log = TRUE) 20.308 21.8925 23.29244 22.6335 23.697 57.334
            stats::dnorm(muKs2, muKs, log = TRUE) 15.832 17.8330 18.98043 18.4695 19.237 52.852
 



 microbenchmark( 
    {mvnfast::dmvn(muKK, rep(0, nKs), covKs, log=TRUE)}, 
    {mvnfast::dmvn(muKK, rep(0, nKs), covKs, log=TRUE,  ncores=8 )},
    {LaplacesDemonCpp::dmvn(muKK, rep(0, nKs), covKs, log=TRUE)},
    {mvtnorm::dmvnorm(muKK, rep(0, nKs), covKs, log=TRUE)},
    {FastGP::rcpp_log_dmvnorm( covKs, rep(0, nKs), muKK, istoep=FALSE )},
    times=1000 )

                                                                      expr      min       lq      mean    median        uq
                {     mvnfast::dmvn(muKK, rep(0, nKs), covKs, log = TRUE) }  191.413  247.684  421.1554  278.4690  302.7620
    {     mvnfast::dmvn(muKK, rep(0, nKs), covKs, log = TRUE, ncores = 8) }  261.157 3134.657 3742.4167 3521.1960 4476.4645
       {     LaplacesDemonCpp::dmvn(muKK, rep(0, nKs), covKs, log = TRUE) }  605.565  665.580 1146.0235  815.9505  890.2885
             {     mvtnorm::dmvnorm(muKK, rep(0, nKs), covKs, log = TRUE) }  591.283  654.936 1119.8288  807.5500  873.9075
 {     FastGP::rcpp_log_dmvnorm(covKs, rep(0, nKs), muKK, istoep = FALSE) } 4840.495 4906.764 5515.3274 5471.8495 6097.7230
        max neval  cld
   3349.043  1000 a   
   6747.828  1000   c 
 115850.571  1000  b  
 115645.814  1000  b  


 nsamp = 1
 microbenchmark( 
    {mvnfast::rmvn(nsamp, rep(0, nKs),  sigma[2]*sigma[2]*exp(-phi*dKK)^kappa )},
    {LaplacesDemon::rmvn(nsamp, rep(0,nKs), sigma[2]*sigma[2]*exp(-phi*dKK)^kappa)},
    {LaplacesDemonCpp::rmvn(nsamp, rep(0,nKs), sigma[2]*sigma[2]*exp(-phi*dKK)^kappa)},
    {mvtnorm::rmvnorm(nsamp, rep(0,nKs), sigma[2]*sigma[2]*exp(-phi*dKK)^kappa, method="chol")},
    {FastGP::rcpp_rmvnorm(nsamp, sigma[2]*sigma[2]*exp(-phi*dKK)^kappa, rep(0,nKs))},
    {FastGP::rcpp_rmvnorm_stable(nsamp, sigma[2]*sigma[2]*exp(-phi*dKK)^kappa, rep(0,nKs))},
    times=1000 )

                                                                                                               expr      min
                     {     mvnfast::rmvn(nsamp, rep(0, nKs), sigma[2] * sigma[2] * exp(-phi *          dKK)^kappa) }  838.675
               {     LaplacesDemon::rmvn(nsamp, rep(0, nKs), sigma[2] * sigma[2] *          exp(-phi * dKK)^kappa) } 6021.432
            {     LaplacesDemonCpp::rmvn(nsamp, rep(0, nKs), sigma[2] * sigma[2] *          exp(-phi * dKK)^kappa) } 2309.990
 {     mvtnorm::rmvnorm(nsamp, rep(0, nKs), sigma[2] * sigma[2] *          exp(-phi * dKK)^kappa, method = "chol") } 1482.535
              {     FastGP::rcpp_rmvnorm(nsamp, sigma[2] * sigma[2] * exp(-phi *          dKK)^kappa, rep(0, nKs)) }  992.656
       {     FastGP::rcpp_rmvnorm_stable(nsamp, sigma[2] * sigma[2] *          exp(-phi * dKK)^kappa, rep(0, nKs)) } 1536.675
       lq     mean   median       uq        max neval  cld
  886.297 1206.807  911.488 1337.527 119179.601  1000 a   
 6193.278 7306.778 6873.217 7745.664 126097.590  1000    d
 2488.021 2768.283 2554.061 3158.878   5174.162  1000   c 
 1564.434 2106.521 1622.486 2291.230 118491.844  1000  b  
 1014.429 1256.690 1036.209 1564.081   3984.498  1000 a   
 1642.857 1991.933 1680.161 2426.348   4620.555  1000  b  
> 

Note:: LaplacesDemonCpp::rmvn is now a copy of mvnfast::rmvn
