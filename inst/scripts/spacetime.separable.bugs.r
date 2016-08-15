
## Bug/Jags model to do spatial-time predition -- separable
## based upon methods of: Viana, M., Jackson, A. L., Graham, N. and Parnell, A. C. 2012. 
## Disentangling spatio-temporal processes in a hierarchical system: a case study in fisheries discards. –
## Ecography 35:
## and Banerjee et al 2008 J. R. Statist. Soc. B (2008) 70, Part 4, pp. 825–848 Gaussian predictive process models for large spatial data sets ( sudipto.bol.ucla.edu/ResearchPapers/BGFS.pdf )

require(sp)
require(rjags)

data(meuse)
data(meuse.grid)

nKs = nrow( meuse  )  # knots 
nPs = nrow( meuse.grid )

xrange = range (c(meuse$x, meuse.grid$x))
yrange = range (c(meuse$y, meuse.grid$y))

dx = diff(xrange)
dy = diff(yrange)
dd = max( dx, dy )

coordsK = data.frame( plon=(meuse$y-yrange[1])/dd, plat=(meuse$x-xrange[1])/dd) 
coordsP = data.frame( plon=(meuse.grid$y-yrange[1])/dd, plat=(meuse.grid$x-xrange[1])/dd) 

# design matrix
X = as.matrix(data.frame( intercept=rep(1, nKs), dist=meuse$dist )) 
XP = as.matrix( data.frame( intercept=rep(1, nPs), dist=meuse.grid$dist ))

y = log(meuse$zinc) 
yP = rep( NA, nPs ) 

dKK = as.matrix(dist( coordsK, diag=TRUE, upper=TRUE)) # distance matrix between knots
dPK = matrix(0, nPs, nKs ) # distance from knot to prediction locations
for (i in 1:nPs) {
  dPK[i,] = sqrt(( coordsK$plon - coordsP$plon[i] )^2 + (coordsK$plat - coordsP$plat[i] )^2)
}

Data = list( nKs=nKs, nPs=nPs, dKK=dKK, dPK=dPK, y=y )

jagsmodel = paste0("
model{
  for(i in 1:nKs){
    y[i] ~ dnorm(mu[i], prec)
    mu[i] = beta0 + errorSpatial[i]
  }
  prec = 1.0/ (tausq + sigmasq.s )
  # Spatial Predictive process
  invcovKs = inverse(covKs)
  errorSpatialK ~ dmnorm( muKs, invcovKs)
  for(i in 1:nKs) {
    muKs[i] = 0
    covKs[i,i] = sigmasq.s
    for(j in 1:(i-1)) {
      covKs[i,j] = sigmasq.s * exp(-(dKK[i,j]/phi.s))
      covKs[j,i] = covKs[i,j]
    } 
  }

  # Interpolate spatial PP back on to original sites
  for(i in 1:nPs) {
    for(j in 1:nKs) {
      covPKs[i,j] = sigmasq.s * exp(-(dPK[i,j]/phi.s))
    } 
  }
  errorSpatial = covPKs %*% invcovKs %*% errorSpatialK
  # errorSpatial = errorSpatialK ## if not using the 'predictive process', this would be used
  
  # Prior distributions
  tausq = 1/tausq.inv
  tausq.inv ~ dgamma(0.1,0.1)
  sigmasq.s = 1/sigmasq.s.inv
  sigmasq.s.inv ~ dgamma(2,1)
  phi.s ~ dgamma(1,0.1)
  beta0 ~ dnorm(0,0.0001)
} 

")

fn = file.path("~/tmp/kriging.bugs")
cat( jagsmodel, file=fn )

fit = jagsUI::jags(data=Data, 
       parameters.to.save=c("phi.s", "sigmasq.s", "tausq"),
       model.file=fn,
       n.iter=1000,
       n.chains=5,
       n.burnin=100,
       n.thin=5,
       parallel=FALSE,
       DIC=FALSE)




## module glm loaded

## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
##    Graph Size: 10216
## 
## Initializing model
This is how the posteriors of ￼ and global.mu look like:

ggs_traceplot(ggs(as.mcmc(fit)))






# -------------------------
# space and time model and some random effects. ..
jagsmodel = paste0("
model{

  for(i in 1:nPs){
    y[i] ∼ dnorm(mu[i], prec[i])
    mu[i] = beta0 + beta1*t[i] + beta2*pow(t[i],2) + u[k[i]] + v[q[i]] + errorSpatial[i] + errorTemporal[i]
    prec[i] = 1/ (tausq + sigmasq.s + sigmasq.t )
  }


  # Spatial Predictive process
  invcovKs = inverse(covKs)
  errorSpatialK ∼ dmnorm(muKs, invcovKs)
  for(i in 1:nKs) {
    muKs[i] = 0
    covKs[i,i] = sigmasq.s
    for(j in 1:(i-1)) {
      covKs[i,j] = sigmasq.s * exp(-(dKK[i,j]/phi.s))
      covKs[j,i] = covKs[i,j]
    } 
  }

  # Interpolate spatial PP back on to original sites
  for(i in 1:nPs) {
    for(j in 1:nKs) {
      covPKs[i,j] = sigmasq.s * exp(-(d.s_k[i,j]/phi.s))
    } 
  }
  errorSpatial = covPKs %*% invcovKs %*% errorSpatialK
  
  # Temporal Predictive process
  errorTemporalK ∼ dmnorm(muKt, invcovKt)
  invcovKt = inverse(covKt)
  for(i in 1:nKt) {
    muKt[i] = 0
    covKt[i,i] = sigmasq.t
    for(j in 1:(i-1)) {
      covKt[i,j] = sigmasq.t * exp(-(d.t_k_k[i,j]/phi.t))
      covKt[j,i] = covKt[i,j]
    } 
  }

  # Interpolate temporal PP back on to original sites
  for(i in 1:nPt) {
    for(j in 1:nKt) {
      covPKt[i,j] = sigmasq.t*exp(-(d.t_k[i,j]/phi.t))   
    } 
  }

  errorTemporal = covPKt %*% invcovKt %*% errorTemporalK

  # Prior distributions
  tausq = 1/tausq.inv
  tausq.inv ∼ dgamma(0.1,0.1)
  sigmasq.s = 1/sigmasq.s.inv
  sigmasq.s.inv ∼ dgamma(2,1)
  phi.s ∼ dgamma(1,0.1)
  beta0 ∼ dnorm(0,0.0001)
  beta1 ∼ dnorm(0,0.0001)
  beta2 ∼ dnorm(0,0.0001)
  sigmasq.t = 1/sigmasq.t.inv
  sigmasq.t.inv ∼ dgamma(2,1)
  phi.t ∼ dunif(7.5,0.5)
  
  # Prior distributions for random effects
  for(i in 1:Nk){ u[i] ∼ dnorm(0,tau.u) }
  tau.u= 1/(sigma.u*sigma.u)
  sigma.u ∼ dunif(0,1000)
  for(i in 1:Nq){ v[i] ∼ dnorm(0,tau.v) }
  tau.v= 1/(sigma.v*sigma.v)
  sigma.v ∼ dunif(0,1000)

} #end model

")



----
 #Random year effects -- AR(1)
      years[1:n_years] ~ dmnorm(zeros_y[],Y.tau[,])
      for(i in 1:n_years){
        for(j in 1:n_years){
          Y.covar[i,j] <- var_year*ifelse(abs(j-i)<=1,pow(rho_year,abs(j-i)),0)
        }
      }
      Y.tau[1:n_years,1:n_years] <- inverse(Y.covar[1:n_years,1:n_years])
      var_year ~ dgamma(1,0.1)
      rho_year ~ dunif(-0.5,0.5)

