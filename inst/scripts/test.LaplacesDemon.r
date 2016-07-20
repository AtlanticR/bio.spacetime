
# ---------------------------------------------------
# basic Universal Kriging
library(sp)
library(gstat)


data(meuse)
coordinates(meuse) = ~x+y

# local universal kriging
gmeuse = gstat(id = "log_zinc", formula = log(zinc)~sqrt(dist), data = meuse)

# variogram of residuals
vmeuse.res = fit.variogram(variogram(gmeuse), vgm(1, "Exp", 300, 1))

# prediction from local neighbourhoods within radius of 170 m or at least 10 points
gmeuse = gstat(id = "log_zinc", formula = log(zinc)~sqrt(dist),
data = meuse, maxdist=170, nmin=10, force=TRUE, model=vmeuse.res)

data(meuse.grid)
gridded(meuse.grid) = ~x+y
predmeuse = predict(gmeuse, meuse.grid)
spplot(predmeuse)



# ---------------------------------------------------
# Universal Kriging with prediction via LaplacesDemon
require(sp)
require(LaplacesDemon)
# require(LaplacesDemonCpp)

data(meuse)
data(meuse.grid)


# Sigma  ## spatial correlation 
# zeta   ## spatial random noise at each location
# beta   ## coeff of linear model 
# X ## linear model (matrix) covariates 
# mu ## mean effect size at each location ( without spatial effect)
# y  ## observations

K = 75 # subset to use to est Sigma

N = nrow( meuse )

xrange = range (c(meuse$x, meuse.grid$x))
yrange = range (c(meuse$y, meuse.grid$y))

dx = diff(xrange)
dy = diff(yrange)
dd = max( dx, dy )

coords = data.frame( plon=(meuse$y-yrange[1])/dd, plat=(meuse$x-xrange[1])/dd) 

coords.new = data.frame( plon=(meuse.grid$y-yrange[1])/dd, plat=(meuse.grid$x-xrange[1])/dd )
  n.new = nrow(coords.new)
  # n.new=1:10
  new = 1:n.new
  coords.new = coords.new[new,]

X = as.matrix(data.frame( intercept=rep(1, N), dist=sqrt( meuse$dist )) )
y = log(meuse$zinc) 

Xnew = as.matrix( data.frame( intercept=rep(1, nrow(meuse.grid)), dist=sqrt( meuse.grid$dist ))) [new,] 
ynew = rep( NA, nrow(meuse.grid) ) [new]

D.N = as.matrix(dist( coords, diag=TRUE, upper=TRUE)) # full distance matrix
D.K = as.matrix(dist( coords[ sample( N, K),], diag=TRUE, upper=TRUE)) # distance matrix

D.new = matrix(0, n.new, K )
for (i in 1:n.new) {
  D.new[i,] = sqrt(( coords$plon[1:K] - coords.new$plon[i] )^2 + (coords$plat[1:K] - coords.new$plat[i] )^2)
}

D.P = matrix(0, N-K, K)
for (i in (K+1):N) {
  D.P[K+1-i,] = sqrt((coords$plon[1:K] - coords$plon[i])^2 + (coords$plat[1:K] - coords$plat[i])^2)
}

mon.names = c("LP",paste0("ynew",new))
parm.names = as.parm.names(list(zeta=rep(0,K), beta=rep(0,2), sigma=rep(0,2), phi=0))
pos.zeta = grep("zeta", parm.names)
pos.beta = grep("beta", parm.names)
pos.sigma = grep("sigma", parm.names)
pos.phi = grep("phi", parm.names)


PGF = function(Data) {
  beta = rnorm(2)
  sigma = runif(2,0.1,10)
  phi = runif(1,1,5)
  kappa = 1
  zeta = mvtnorm::rmvnorm(1, rep(0,Data$K), sigma[2]*sigma[2]*exp(-phi*Data$D.K)^kappa, method="chol")
  return(c(zeta, beta, sigma, phi))
}


Data = list(D.N=D.N, D.K=D.K, D.new=D.new, D.P=D.P, K=K, N=N, PGF=PGF, X=X,
  Xnew=Xnew, plat=coords$plat, plon=coords$plon,
  mon.names=mon.names, parm.names=parm.names, pos.zeta=pos.zeta,
  pos.beta=pos.beta, pos.sigma=pos.sigma, pos.phi=pos.phi, y=y)


Model = function(parm, Data){
  ### Parameters
  beta = parm[Data$pos.beta]
  zeta = parm[Data$pos.zeta]
  kappa = 1
  parm[Data$pos.sigma] = sigma = interval(parm[Data$pos.sigma], 1, Inf)
  parm[Data$pos.phi] = phi = interval(parm[Data$pos.phi], 1, 5)
  Sigma = sigma[2]*sigma[2] * exp(-phi * Data$D.K)^kappa
  
  ### Log-Prior
  beta.prior = sum(dnormv(beta, 0, 1000, log=TRUE))
  zeta.prior = mvtnorm::dmvnorm(zeta, rep(0, Data$K), Sigma, log=TRUE)
  sigma.prior = sum(dhalfcauchy(sigma - 1, 25, log=TRUE))
  phi.prior = dunif(phi, 1, 5, log=TRUE)
  
  ### Interpolation
  rho = exp(-phi * Data$D.new)^kappa   ## spatial correlation
  ynew = rnorm(nrow(Data$Xnew), tcrossprod(Data$Xnew, t(beta)) + rowSums(rho / rowSums(rho) * zeta), sigma[1]) 

  ### Log-Likelihood
  mu = tcrossprod(Data$X, t(beta))
  mu[1:Data$K] = mu[1:Data$K] + zeta

  lambda = exp(-phi * Data$D.P)^kappa
  ii = (Data$K+1):Data$N
  mu[ii] = mu[ii] + rowSums(lambda / rowSums(lambda) * matrix(zeta, Data$N - Data$K, Data$K, byrow=TRUE))
  LL = sum(dnorm(Data$y, mu, sigma[1], log=TRUE))
  ### Log-Posterior
  LP = LL + beta.prior + zeta.prior + sigma.prior + phi.prior
  Modelout = list(LP=LP, Dev=-2*LL, Monitor=c(LP,ynew), yhat=rnorm(length(mu), mu, sigma[1]), parm=parm)
  return(Modelout)
}

# Initial.Values = c(rep(0,K), c(mean(y), 0), rep(1,2), 3)

f = LaplaceApproximation(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000, Method="SPG" ) # fast inital solution
f = LaplaceApproximation(Model, Data=Data, parm=as.initial.values(f), Iterations=1000, Method="LBFGS" ) # refine it
f = VariationalBayes(Model, Data=Data, parm=as.initial.values(f), Iterations=500, Samples=20, CPUs=5, Covar=f$Covar ) # refine it again


f = LaplacesDemon(Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=1000, Status=100, Thinning=1, Covar=f$Covar )
f = VariationalBayes(Model, Data=Data, parm=as.initial.values(f), Iterations=100, Samples=10, CPUs=5, Covar=f$Covar )
f = IterativeQuadrature(Model, Data=Data, parm=as.initial.values(f), Iterations=10, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL), Covar=f$Covar )


inew = grep( "ynew", rownames( f$Summary2 ) )
m = f$Summary2[inew,]
preds = coords.new 
preds$prediction = f$Summary2[inew, 1]




---
# space-time separable
data(demontexas)
Y = as.matrix(demontexas[1:20,c(18:30)])
X = cbind(1,as.matrix(demontexas[1:20,c(1,4)])) #Static predictors
plat = demontexas[1:20,2]
plon = demontexas[1:20,3]
S = nrow(Y) #Number of sites, or points in space
T = ncol(Y) #Number of time-periods
K = ncol(X) #Number of columns in design matrix X including the intercept
D.S = as.matrix(dist(cbind(plon,plat), diag=TRUE, upper=TRUE))
D.T = as.matrix(dist(cbind(c(1:T),c(1:T)), diag=TRUE, upper=TRUE))
mon.names = "LP"
parm.names = as.parm.names(list(zeta=rep(0,S), theta=rep(0,T),
beta=rep(0,K), phi=rep(0,2), sigma=rep(0,3)))
pos.zeta = grep("zeta", parm.names)
pos.theta = grep("theta", parm.names)
pos.beta = grep("beta", parm.names)
pos.phi = grep("phi", parm.names)
pos.sigma = grep("sigma", parm.names)

PGF = function(Data) {
  beta = rnorm(Data$K, c(mean(Data$Y),rep(0,Data$K-1)), 1)
  phi = runif(2,1,5)
  sigma = runif(3)
  kappa = 1
  lambda = 1
  Sigma.S = sigma[2]^2 * exp(-phi[1] * Data$D.S)^kappa
  Sigma.T = sigma[3]^2 * exp(-phi[2] * Data$D.T)^lambda
  zeta = as.vector(mvtnorm::rmvnorm(1, rep(0,Data$S), Sigma.S, method="chol" ))
  theta = as.vector(mvtnorm::rmvnorm(1, rep(0,Data$T), Sigma.T, method="chol" ))
  return(c(zeta, theta, beta, phi, sigma))
}

Data = list(D.S=D.S, D.T=D.T, K=K, PGF=PGF, S=S, T=T, X=X, Y=Y,
  plat=plat, plon=plon, mon.names=mon.names,
  parm.names=parm.names, pos.zeta=pos.zeta, pos.theta=pos.theta,
  pos.beta=pos.beta, pos.phi=pos.phi, pos.sigma=pos.sigma)

Model = function(parm, Data) {
  ### Hyperparameters
  zeta.mu = rep(0,Data$S)
  theta.mu = rep(0,Data$T)
  ### Parameters
  beta = parm[Data$pos.beta]
  zeta = parm[Data$pos.zeta]
  theta = parm[Data$pos.theta]
  kappa = 1; lambda = 1
  sigma = interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] = sigma
  parm[Data$pos.phi] = phi = interval(parm[Data$pos.phi], 1, 5)
  Sigma.S = sigma[2]^2 * exp(-phi[1] * Data$D.S)^kappa
  Sigma.T = sigma[3]^2 * exp(-phi[2] * Data$D.T)^lambda
  ### Log-Prior
  beta.prior = sum(dnormv(beta, 0, 1000, log=TRUE))
  zeta.prior = mvtnorm::dmvnorm(zeta, zeta.mu, Sigma.S, log=TRUE)
  theta.prior = mvtnorm::dmvnorm(theta, theta.mu, Sigma.T, log=TRUE)
  sigma.prior = sum(dhalfcauchy(sigma, 25, log=TRUE))
  phi.prior = sum(dunif(phi, 1, 5, log=TRUE))
  ### Log-Likelihood
  Theta = matrix(theta, Data$S, Data$T, byrow=TRUE)
  mu = as.vector(tcrossprod(Data$X, t(beta))) + zeta + Theta
  LL = sum(dnorm(Data$Y, mu, sigma[1], log=TRUE))
  ### Log-Posterior
  LP = LL + beta.prior + zeta.prior + theta.prior + sigma.prior +
  phi.prior
  Modelout = list(LP=LP, Dev=-2*LL, Monitor=LP,
  yhat=rnorm(prod(dim(mu)), mu, sigma[1]), parm=parm)
  return(Modelout)
}

Initial.Values = c(rep(0,S), rep(0,T), rep(0,2), rep(1,2), rep(1,3))


f = LaplaceApproximation(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000 )
Initial.Values = as.initial.values(f)
f = LaplacesDemon(Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=1000, Status=100, Thinning=1 )
f = VariationalBayes(Model, Data=Data, parm=as.initial.values(f), Iterations=10000, Samples=1000, CPUs=5 )
f = IterativeQuadrature(Model, Data=Data, parm=as.initial.values(f), Iterations=1000, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )



---

# space time non-separable

data(demontexas)
Y = as.matrix(demontexas[1:10,c(18:30)])
X = cbind(1,as.matrix(demontexas[1:10,c(1,4)])) #Static predictors
plat = demontexas[1:10,2]
plon = demontexas[1:10,3]
S = nrow(Y) #Number of sites, or points in space
T = ncol(Y) #Number of time-periods
K = ncol(X) #Number of columns in design matrix X including the intercept
D.S = as.matrix(dist(cbind(rep(plon,T),rep(plat,T)), diag=TRUE,
upper=TRUE))
D.T = as.matrix(dist(cbind(rep(1:T,each=S),rep(1:T,each=S)), diag=TRUE,
upper=TRUE))
mon.names = "LP"
parm.names = as.parm.names(list(Xi=matrix(0,S,T), beta=rep(0,K),
phi=rep(0,2), sigma=rep(0,2), psi=0))
pos.Xi = grep("Xi", parm.names)
pos.beta = grep("beta", parm.names)
pos.phi = grep("phi", parm.names)
pos.sigma = grep("sigma", parm.names)
pos.psi = grep("psi", parm.names)
PGF = function(Data) {
beta = rnorm(Data$K, c(mean(Data$Y),rep(0,Data$K-1)), 1)
phi = runif(2,1,5)
sigma = runif(2)
psi = runif(1)
kappa = 1
lambda = 1
Sigma = sigma[2]*sigma[2] * exp(-(Data$D.S / phi[1])^kappa - (Data$D.T
/ phi[2])^lambda - psi*(Data$D.S / phi[1])^kappa * (Data$D.T / phi[2])^lambda)
Xi = as.vector(rmvn(1, rep(0,Data$S*Data$T), Sigma))
return(c(Xi, beta, phi, sigma, psi))
}
Data = list(D.S=D.S, D.T=D.T, K=K, PGF=PGF, S=S, T=T, X=X, Y=Y,
plat=plat, plon=plon, mon.names=mon.names,
parm.names=parm.names, pos.Xi=pos.Xi, pos.beta=pos.beta,
pos.phi=pos.phi, pos.sigma=pos.sigma, pos.psi=pos.psi)

Model = function(parm, Data) {
  ### Hyperparameters
  Xi.mu = rep(0,Data$S*Data$T)
  ### Parameters
  beta = parm[Data$pos.beta]
  Xi = parm[Data$pos.Xi]
  kappa = 1; lambda = 1
  sigma = interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] = sigma
  parm[Data$pos.phi] = phi = interval(parm[Data$pos.phi], 1, 5)
  parm[Data$pos.psi] = psi = interval(parm[Data$pos.psi], 1e-100, Inf)
  Sigma = sigma[2]*sigma[2] * exp(-(Data$D.S / phi[1])^kappa -
  (Data$D.T / phi[2])^lambda -
  psi*(Data$D.S / phi[1])^kappa * (Data$D.T / phi[2])^lambda)
  ### Log-Prior
  beta.prior = sum(dnormv(beta, 0, 1000, log=TRUE))
  Xi.prior = mvtnorm::dmvnorm(Xi, Xi.mu, Sigma, log=TRUE)
  sigma.prior = sum(dhalfcauchy(sigma, 25, log=TRUE))
  phi.prior = sum(dunif(phi, 1, 5, log=TRUE))
  psi.prior = dhalfcauchy(psi, 25, log=TRUE)
  ### Log-Likelihood
  Xi = matrix(Xi, Data$S, Data$T)
  mu = as.vector(tcrossprod(Data$X, t(beta))) + Xi
  LL = sum(dnorm(Data$Y, mu, sigma[1], log=TRUE))
  ### Log-Posterior
  LP = LL + beta.prior + Xi.prior + sigma.prior + phi.prior + psi.prior
  Modelout = list(LP=LP, Dev=-2*LL, Monitor=LP,
  yhat=rnorm(prod(dim(mu)), mu, sigma[1]), parm=parm)
  return(Modelout)
}

Initial.Values = c(rep(0,S*T), c(mean(Y),rep(0,K-1)), rep(1,2), rep(1,2), 1)


f = LaplaceApproximation(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000 )
f = LaplacesDemon(Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=1000, Status=100, Thinning=1 )
f = VariationalBayes(Model, Data=Data, parm=as.initial.values(f), Iterations=10000, Samples=1000, CPUs=5 )
f = IterativeQuadrature(Model, Data=Data, parm=as.initial.values(f), Iterations=1000, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )




----


# spacetime dynamic
data(demontexas)
Y = as.matrix(demontexas[1:20,c(18:30)])
X = cbind(1,as.matrix(demontexas[1:20,c(1,4)])) #Static predictors
plat = demontexas[1:20,2]
plon = demontexas[1:20,3]
D = as.matrix(dist(cbind(plon,plat), diag=TRUE, upper=TRUE))
S = nrow(Y) #Number of sites, or points in space
T = ncol(Y) #Number of time-periods
K = ncol(X) #Number of columns in design matrix X including the intercept
mon.names = "LP"
parm.names = as.parm.names(list(zeta=rep(0,S), beta=matrix(0,K,T),
phi=rep(0,T), kappa=rep(0,T), lambda=rep(0,T), sigma=rep(0,4),
tau=rep(0,K)))
pos.zeta = grep("zeta", parm.names)
pos.beta = grep("beta", parm.names)
pos.phi = grep("phi", parm.names)
pos.kappa = grep("kappa", parm.names)
pos.lambda = grep("lambda", parm.names)
pos.sigma = grep("sigma", parm.names)

pos.tau = grep("tau", parm.names)
PGF = function(Data) {
  beta = rnorm(Data$K*Data$T, rbind(mean(Data$Y),
  matrix(0, Data$K-1, Data$T)), 1)
  phi = rhalfnorm(Data$T, 1)
  kappa = rhalfnorm(Data$T, 1)
  lambda = rhalfnorm(Data$T, 1)
  Sigma = lambda[1]*lambda[1]*exp(-phi[1]*Data$D)^kappa[1]
  zeta = as.vector(mvtnorm::rmvnorm(1, rep(0,Data$S), Sigma, method="chol" ))
  sigma = runif(4)
  tau = runif(Data$K)
  return(c(zeta, beta, phi, kappa, lambda, sigma, tau))
}

Data = list(D=D, K=K, PGF=PGF, S=S, T=T, X=X, Y=Y, plat=plat,
  plon=plon, mon.names=mon.names, parm.names=parm.names,
  pos.zeta=pos.zeta, pos.beta=pos.beta, pos.phi=pos.phi,
  pos.kappa=pos.kappa, pos.lambda=pos.lambda, pos.sigma=pos.sigma,
  pos.tau=pos.tau)

Model = function(parm, Data) {
  ### Parameters
  beta = matrix(parm[Data$pos.beta], Data$K, Data$T)
  zeta = parm[Data$pos.zeta]
  parm[Data$pos.phi] = phi = interval(parm[Data$pos.phi], 1e-100, Inf)
  kappa = interval(parm[Data$pos.kappa], 1e-100, Inf)
  parm[Data$pos.kappa] = kappa
  lambda = interval(parm[Data$pos.lambda], 1e-100, Inf)
  parm[Data$pos.lambda] = lambda
  sigma = interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] = sigma
  parm[Data$pos.tau] = tau = interval(parm[Data$pos.tau], 1e-100, Inf)
  Sigma = array(0, dim=c(Data$S, Data$S, Data$T))
  for (t in 1:Data$T) {
  Sigma[ , ,t] = lambda[t]^2 * exp(-phi[t] * Data$D)^kappa[t]}
  ### Log-Prior
  beta.prior = sum(dnormv(beta[,1], 0, 1000, log=TRUE),
  dnorm(beta[,-1], beta[,-Data$T], matrix(tau, Data$K,
  Data$T-1), log=TRUE))
  zeta.prior = mvtnorm::dmvnorm(zeta, rep(0,Data$S), Sigma[ , , 1], log=TRUE)
  phi.prior = sum(dhalfnorm(phi[1], sqrt(1000), log=TRUE),
  dtrunc(phi[-1], "norm", a=0, b=Inf, mean=phi[-Data$T],
  sd=sigma[2], log=TRUE))
  kappa.prior = sum(dhalfnorm(kappa[1], sqrt(1000), log=TRUE),
  dtrunc(kappa[-1], "norm", a=0, b=Inf, mean=kappa[-Data$T],
  sd=sigma[3], log=TRUE))
  lambda.prior = sum(dhalfnorm(lambda[1], sqrt(1000), log=TRUE),
  dtrunc(lambda[-1], "norm", a=0, b=Inf, mean=lambda[-Data$T],
  sd=sigma[4], log=TRUE))
  sigma.prior = sum(dhalfcauchy(sigma, 25, log=TRUE))
  tau.prior = sum(dhalfcauchy(tau, 25, log=TRUE))
  ### Log-Likelihood
  mu = tcrossprod(Data$X, t(beta))
  Theta = matrix(zeta, Data$S, Data$T)
  for (t in 2:Data$T) {
  for (s in 1:Data$S) {
  Theta[s,t] = sum(Sigma[,s,t] / sum(Sigma[,s,t]) * Theta[,t-1])}}
  mu = mu + Theta
  LL = sum(dnorm(Data$Y, mu, sigma[1], log=TRUE))
  ### Log-Posterior
  LP = LL + beta.prior + zeta.prior + sum(phi.prior) +
  sum(kappa.prior) + sum(lambda.prior) + sigma.prior + tau.prior
  Modelout = list(LP=LP, Dev=-2*LL, Monitor=LP,
  yhat=rnorm(prod(dim(mu)), mu, sigma[1]), parm=parm)
  return(Modelout)
}

Initial.Values = c(rep(0,S), rep(c(mean(Y),rep(0,K-1)),T), rep(1,T),
rep(1,T), rep(1,T), rep(1,4), rep(1,K))


f = LaplaceApproximation(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000 )

f = LaplacesDemon(Model, Data=Data, Initial.Values=as.initial.values(f), Iterations=1000, Status=100, Thinning=1 )
f = VariationalBayes(Model, Data=Data, parm=as.initial.values(f), Iterations=10000, Samples=1000, CPUs=5 )
f = IterativeQuadrature(Model, Data=Data, parm=as.initial.values(f), Iterations=1000, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )

