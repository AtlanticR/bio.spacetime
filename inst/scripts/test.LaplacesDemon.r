
# basic Univ Kriging

require(LaplacesDemon)
require(LaplacesDemonCpp)


N <- 20
longitude <- runif(N+1,0,100)
latitude <- runif(N+1,0,100)
D <- as.matrix(dist(cbind(longitude,latitude), diag=TRUE, upper=TRUE))
Sigma <- 10000 * exp(-1.5 * D)
zeta <- colMeans(rmvn(1000, rep(0,N+1), Sigma))
beta <- c(50,2)
X <- matrix(runif((N+1)*2,-2,2),(N+1),2); X[,1] <- 1
mu <- as.vector(tcrossprod(X, t(beta)))
y <- mu + zeta
longitude.new <- longitude[N+1]; latitude.new <- latitude[N+1]
Xnew <- X[N+1,]; ynew <- y[N+1]
longitude <- longitude[1:N]; latitude <- latitude[1:N]
X <- X[1:N,]; y <- y[1:N]
D <- as.matrix(dist(cbind(longitude,latitude), diag=TRUE, upper=TRUE))
D.new <- sqrt((longitude - longitude.new)^2 + (latitude - latitude.new)^2)
mon.names <- c("LP","ynew")
parm.names <- as.parm.names(list(zeta=rep(0,N), beta=rep(0,2),
sigma=rep(0,2), phi=0))
pos.zeta <- grep("zeta", parm.names)
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.phi <- grep("phi", parm.names)

PGF <- function(Data) {
  beta <- rnorm(2)
  sigma <- runif(2,0.1,10)
  phi <- runif(1,1,5)
  kappa <- 1
  zeta <- rmvn(1, rep(0,Data$N),sigma[2]*sigma[2]*exp(-phi*Data$D)^kappa)
  return(c(zeta, beta, sigma, phi))
}


Data <- list(D=D, D.new=D.new, latitude=latitude, longitude=longitude,
  N=N, PGF=PGF, X=X, Xnew=Xnew, mon.names=mon.names,
  parm.names=parm.names, pos.zeta=pos.zeta, pos.beta=pos.beta,
  pos.sigma=pos.sigma, pos.phi=pos.phi, y=y)


Model <- function(parm, Data){
  ### Parameters
  beta <- parm[Data$pos.beta]
  zeta <- parm[Data$pos.zeta]
  kappa <- 1
  sigma <- interval(parm[Data$pos.sigma], 0.1, 10)
  parm[Data$pos.sigma] <- sigma
  parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], 1, 5)
  Sigma <- sigma[2]*sigma[2] * exp(-phi * Data$D)^kappa
  ### Log-Prior
  beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
  zeta.prior <- dmvn(zeta, rep(0, Data$N), Sigma, log=TRUE)
  sigma.prior <- sum(dhalfcauchy(sigma - 1, 25, log=TRUE))
  phi.prior <- dunif(phi, 1, 5, log=TRUE)
  ### Interpolation
  rho <- exp(-phi * Data$D.new)^kappa
  ynew <- rnorm(1, sum(beta * Data$Xnew) + sum(rho / sum(rho) * zeta),
  sigma[1])
  ### Log-Likelihood
  mu <- tcrossprod(Data$X, t(beta)) + zeta
  LL <- sum(dnorm(Data$y, mu, sigma[1], log=TRUE))
  ### Log-Posterior
  LP <- LL + beta.prior + zeta.prior + sigma.prior + phi.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,ynew),
  yhat=rnorm(length(mu), mu, sigma[1]), parm=parm)
  return(Modelout)
}

Initial.Values <- c(rep(0,N), rep(0,2), rep(1,2), 1)

Fit = LaplaceApproximation(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000 )
Initial.Values <- as.initial.values(Fit)
Fit = LaplacesDemon(Model, Data=Data, Initial.Values=Initial.Values, Iterations=1000, Status=100, Thinning=1 )
Fit = VariationalBayes(Model, Data=Data, parm=Initial.Values, Iterations=10000, Samples=1000, CPUs=5 )
Fit = IterativeQuadrature(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )





-----
# UK with prediction

N <- 100
K <- 30 #Number of knots
longitude <- runif(N+1,0,100)
latitude <- runif(N+1,0,100)
D <- as.matrix(dist(cbind(longitude,latitude), diag=TRUE, upper=TRUE))
Sigma <- 10000 * exp(-1.5 * D)
zeta <- colMeans(rmvn(1000, rep(0,N+1), Sigma))
beta <- c(50,2)
X <- matrix(runif((N+1)*2,-2,2),(N+1),2); X[,1] <- 1
mu <- as.vector(tcrossprod(X, t(beta)))
y <- mu + zeta
longitude.new <- longitude[N+1]; latitude.new <- latitude[N+1]
Xnew <- X[N+1,]; ynew <- y[N+1]
longitude <- longitude[1:N]; latitude <- latitude[1:N]
X <- X[1:N,]; y <- y[1:N]
D <- as.matrix(dist(cbind(longitude[1:K],latitude[1:K]), diag=TRUE,
upper=TRUE))
D.P <- matrix(0, N-K, K)

for (i in (K+1):N) {
  D.P[K+1-i,] <- sqrt((longitude[1:K] - longitude[i])^2 + (latitude[1:K] - latitude[i])^2)
}
D.new <- sqrt((longitude[1:K] - longitude.new)^2 +
(latitude[1:K] - latitude.new)^2)
mon.names <- c("LP","ynew")
parm.names <- as.parm.names(list(zeta=rep(0,K), beta=rep(0,2),
sigma=rep(0,2), phi=0))
pos.zeta <- grep("zeta", parm.names)
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.phi <- grep("phi", parm.names)

PGF <- function(Data) {
  beta <- rnorm(2)
  sigma <- runif(2,0.1,10)
  phi <- runif(1,1,5)
  kappa <- 1
  zeta <- rmvn(1, rep(0,Data$K),
  sigma[2]*sigma[2]*exp(-phi*Data$D)^kappa)
  return(c(zeta, beta, sigma, phi))
}

Data <- list(D=D, D.new=D.new, D.P=D.P, K=K, N=N, PGF=PGF, X=X,
  Xnew=Xnew, latitude=latitude, longitude=longitude,
  mon.names=mon.names, parm.names=parm.names, pos.zeta=pos.zeta,
  pos.beta=pos.beta, pos.sigma=pos.sigma, pos.phi=pos.phi, y=y)


Model <- function(parm, Data){
  ### Parameters
  beta <- parm[Data$pos.beta]
  zeta <- parm[Data$pos.zeta]
  kappa <- 1
  sigma <- interval(parm[Data$pos.sigma], 1, Inf)
  parm[Data$pos.sigma] <- sigma
  parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], 1, 5)
  Sigma <- sigma[2]*sigma[2] * exp(-phi * Data$D)^kappa
  ### Log-Prior
  beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
  zeta.prior <- dmvn(zeta, rep(0, Data$K), Sigma, log=TRUE)
  sigma.prior <- sum(dhalfcauchy(sigma - 1, 25, log=TRUE))
  phi.prior <- dunif(phi, 1, 5, log=TRUE)
  ### Interpolation
  rho <- exp(-phi * Data$D.new)^kappa
  ynew <- rnorm(1, sum(beta * Data$Xnew) + sum(rho / sum(rho) * zeta),sigma)

  ### Log-Likelihood
  mu <- tcrossprod(Data$X, t(beta))
  mu[1:Data$K] <- mu[1:Data$K] + zeta
  lambda <- exp(-phi * Data$D.P)^kappa
  mu[(Data$K+1):Data$N] <- mu[(Data$K+1):Data$N] +
  rowSums(lambda / rowSums(lambda) *
  matrix(zeta, Data$N - Data$K, Data$K, byrow=TRUE))
  LL <- sum(dnorm(Data$y, mu, sigma[1], log=TRUE))
  ### Log-Posterior
  LP <- LL + beta.prior + zeta.prior + sigma.prior + phi.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,ynew),
  yhat=rnorm(length(mu), mu, sigma[1]), parm=parm)
  return(Modelout)
}

Initial.Values <- c(rep(0,K), c(mean(y), 0), rep(1,2), 3)

Fit = LaplaceApproximation(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000 )
Initial.Values <- as.initial.values(Fit)
Fit = LaplacesDemon(Model, Data=Data, Initial.Values=Initial.Values, Iterations=1000, Status=100, Thinning=1 )
Fit = VariationalBayes(Model, Data=Data, parm=Initial.Values, Iterations=10000, Samples=1000, CPUs=5 )
Fit = IterativeQuadrature(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )





----

# space-time separable
data(demontexas)
Y <- as.matrix(demontexas[1:20,c(18:30)])
X <- cbind(1,as.matrix(demontexas[1:20,c(1,4)])) #Static predictors
latitude <- demontexas[1:20,2]
longitude <- demontexas[1:20,3]
S <- nrow(Y) #Number of sites, or points in space
T <- ncol(Y) #Number of time-periods
K <- ncol(X) #Number of columns in design matrix X including the intercept
D.S <- as.matrix(dist(cbind(longitude,latitude), diag=TRUE, upper=TRUE))
D.T <- as.matrix(dist(cbind(c(1:T),c(1:T)), diag=TRUE, upper=TRUE))
mon.names <- "LP"
parm.names <- as.parm.names(list(zeta=rep(0,S), theta=rep(0,T),
beta=rep(0,K), phi=rep(0,2), sigma=rep(0,3)))
pos.zeta <- grep("zeta", parm.names)
pos.theta <- grep("theta", parm.names)
pos.beta <- grep("beta", parm.names)
pos.phi <- grep("phi", parm.names)
pos.sigma <- grep("sigma", parm.names)

PGF <- function(Data) {
  beta <- rnorm(Data$K, c(mean(Data$Y),rep(0,Data$K-1)), 1)
  phi <- runif(2,1,5)
  sigma <- runif(3)
  kappa <- 1
  lambda <- 1
  Sigma.S <- sigma[2]^2 * exp(-phi[1] * Data$D.S)^kappa
  Sigma.T <- sigma[3]^2 * exp(-phi[2] * Data$D.T)^lambda
  zeta <- as.vector(rmvn(1, rep(0,Data$S), Sigma.S))
  theta <- as.vector(rmvn(1, rep(0,Data$T), Sigma.T))
  return(c(zeta, theta, beta, phi, sigma))
}

Data <- list(D.S=D.S, D.T=D.T, K=K, PGF=PGF, S=S, T=T, X=X, Y=Y,
  latitude=latitude, longitude=longitude, mon.names=mon.names,
  parm.names=parm.names, pos.zeta=pos.zeta, pos.theta=pos.theta,
  pos.beta=pos.beta, pos.phi=pos.phi, pos.sigma=pos.sigma)

Model <- function(parm, Data) {
  ### Hyperparameters
  zeta.mu <- rep(0,Data$S)
  theta.mu <- rep(0,Data$T)
  ### Parameters
  beta <- parm[Data$pos.beta]
  zeta <- parm[Data$pos.zeta]
  theta <- parm[Data$pos.theta]
  kappa <- 1; lambda <- 1
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma
  parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], 1, 5)
  Sigma.S <- sigma[2]^2 * exp(-phi[1] * Data$D.S)^kappa
  Sigma.T <- sigma[3]^2 * exp(-phi[2] * Data$D.T)^lambda
  ### Log-Prior
  beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
  zeta.prior <- dmvn(zeta, zeta.mu, Sigma.S, log=TRUE)
  theta.prior <- dmvn(theta, theta.mu, Sigma.T, log=TRUE)
  sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
  phi.prior <- sum(dunif(phi, 1, 5, log=TRUE))
  ### Log-Likelihood
  Theta <- matrix(theta, Data$S, Data$T, byrow=TRUE)
  mu <- as.vector(tcrossprod(Data$X, t(beta))) + zeta + Theta
  LL <- sum(dnorm(Data$Y, mu, sigma[1], log=TRUE))
  ### Log-Posterior
  LP <- LL + beta.prior + zeta.prior + theta.prior + sigma.prior +
  phi.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
  yhat=rnorm(prod(dim(mu)), mu, sigma[1]), parm=parm)
  return(Modelout)
}

Initial.Values <- c(rep(0,S), rep(0,T), rep(0,2), rep(1,2), rep(1,3))


Fit = LaplaceApproximation(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000 )
Initial.Values <- as.initial.values(Fit)
Fit = LaplacesDemon(Model, Data=Data, Initial.Values=Initial.Values, Iterations=1000, Status=100, Thinning=1 )
Fit = VariationalBayes(Model, Data=Data, parm=Initial.Values, Iterations=10000, Samples=1000, CPUs=5 )
Fit = IterativeQuadrature(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )



---

# space time non-separable

data(demontexas)
Y <- as.matrix(demontexas[1:10,c(18:30)])
X <- cbind(1,as.matrix(demontexas[1:10,c(1,4)])) #Static predictors
latitude <- demontexas[1:10,2]
longitude <- demontexas[1:10,3]
S <- nrow(Y) #Number of sites, or points in space
T <- ncol(Y) #Number of time-periods
K <- ncol(X) #Number of columns in design matrix X including the intercept
D.S <- as.matrix(dist(cbind(rep(longitude,T),rep(latitude,T)), diag=TRUE,
upper=TRUE))
D.T <- as.matrix(dist(cbind(rep(1:T,each=S),rep(1:T,each=S)), diag=TRUE,
upper=TRUE))
mon.names <- "LP"
parm.names <- as.parm.names(list(Xi=matrix(0,S,T), beta=rep(0,K),
phi=rep(0,2), sigma=rep(0,2), psi=0))
pos.Xi <- grep("Xi", parm.names)
pos.beta <- grep("beta", parm.names)
pos.phi <- grep("phi", parm.names)
pos.sigma <- grep("sigma", parm.names)
pos.psi <- grep("psi", parm.names)
PGF <- function(Data) {
beta <- rnorm(Data$K, c(mean(Data$Y),rep(0,Data$K-1)), 1)
phi <- runif(2,1,5)
sigma <- runif(2)
psi <- runif(1)
kappa <- 1
lambda <- 1
Sigma <- sigma[2]*sigma[2] * exp(-(Data$D.S / phi[1])^kappa - (Data$D.T
/ phi[2])^lambda - psi*(Data$D.S / phi[1])^kappa * (Data$D.T / phi[2])^lambda)
Xi <- as.vector(rmvn(1, rep(0,Data$S*Data$T), Sigma))
return(c(Xi, beta, phi, sigma, psi))
}
Data <- list(D.S=D.S, D.T=D.T, K=K, PGF=PGF, S=S, T=T, X=X, Y=Y,
latitude=latitude, longitude=longitude, mon.names=mon.names,
parm.names=parm.names, pos.Xi=pos.Xi, pos.beta=pos.beta,
pos.phi=pos.phi, pos.sigma=pos.sigma, pos.psi=pos.psi)

Model <- function(parm, Data) {
  ### Hyperparameters
  Xi.mu <- rep(0,Data$S*Data$T)
  ### Parameters
  beta <- parm[Data$pos.beta]
  Xi <- parm[Data$pos.Xi]
  kappa <- 1; lambda <- 1
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma
  parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], 1, 5)
  parm[Data$pos.psi] <- psi <- interval(parm[Data$pos.psi], 1e-100, Inf)
  Sigma <- sigma[2]*sigma[2] * exp(-(Data$D.S / phi[1])^kappa -
  (Data$D.T / phi[2])^lambda -
  psi*(Data$D.S / phi[1])^kappa * (Data$D.T / phi[2])^lambda)
  ### Log-Prior
  beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
  Xi.prior <- dmvn(Xi, Xi.mu, Sigma, log=TRUE)
  sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
  phi.prior <- sum(dunif(phi, 1, 5, log=TRUE))
  psi.prior <- dhalfcauchy(psi, 25, log=TRUE)
  ### Log-Likelihood
  Xi <- matrix(Xi, Data$S, Data$T)
  mu <- as.vector(tcrossprod(Data$X, t(beta))) + Xi
  LL <- sum(dnorm(Data$Y, mu, sigma[1], log=TRUE))
  ### Log-Posterior
  LP <- LL + beta.prior + Xi.prior + sigma.prior + phi.prior + psi.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
  yhat=rnorm(prod(dim(mu)), mu, sigma[1]), parm=parm)
  return(Modelout)
}

Initial.Values <- c(rep(0,S*T), c(mean(Y),rep(0,K-1)), rep(1,2), rep(1,2), 1)


Fit = LaplaceApproximation(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000 )
Initial.Values <- as.initial.values(Fit)
Fit = LaplacesDemon(Model, Data=Data, Initial.Values=Initial.Values, Iterations=1000, Status=100, Thinning=1 )
Fit = VariationalBayes(Model, Data=Data, parm=Initial.Values, Iterations=10000, Samples=1000, CPUs=5 )
Fit = IterativeQuadrature(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )




----


# spacetime dynamic
data(demontexas)
Y <- as.matrix(demontexas[1:20,c(18:30)])
X <- cbind(1,as.matrix(demontexas[1:20,c(1,4)])) #Static predictors
latitude <- demontexas[1:20,2]
longitude <- demontexas[1:20,3]
D <- as.matrix(dist(cbind(longitude,latitude), diag=TRUE, upper=TRUE))
S <- nrow(Y) #Number of sites, or points in space
T <- ncol(Y) #Number of time-periods
K <- ncol(X) #Number of columns in design matrix X including the intercept
mon.names <- "LP"
parm.names <- as.parm.names(list(zeta=rep(0,S), beta=matrix(0,K,T),
phi=rep(0,T), kappa=rep(0,T), lambda=rep(0,T), sigma=rep(0,4),
tau=rep(0,K)))
pos.zeta <- grep("zeta", parm.names)
pos.beta <- grep("beta", parm.names)
pos.phi <- grep("phi", parm.names)
pos.kappa <- grep("kappa", parm.names)
pos.lambda <- grep("lambda", parm.names)
pos.sigma <- grep("sigma", parm.names)

pos.tau <- grep("tau", parm.names)
PGF <- function(Data) {
  beta <- rnorm(Data$K*Data$T, rbind(mean(Data$Y),
  matrix(0, Data$K-1, Data$T)), 1)
  phi <- rhalfnorm(Data$T, 1)
  kappa <- rhalfnorm(Data$T, 1)
  lambda <- rhalfnorm(Data$T, 1)
  Sigma <- lambda[1]*lambda[1]*exp(-phi[1]*Data$D)^kappa[1]
  zeta <- as.vector(rmvn(1, rep(0,Data$S), Sigma))
  sigma <- runif(4)
  tau <- runif(Data$K)
  return(c(zeta, beta, phi, kappa, lambda, sigma, tau))
}

Data <- list(D=D, K=K, PGF=PGF, S=S, T=T, X=X, Y=Y, latitude=latitude,
  longitude=longitude, mon.names=mon.names, parm.names=parm.names,
  pos.zeta=pos.zeta, pos.beta=pos.beta, pos.phi=pos.phi,
  pos.kappa=pos.kappa, pos.lambda=pos.lambda, pos.sigma=pos.sigma,
  pos.tau=pos.tau)

Model <- function(parm, Data) {
  ### Parameters
  beta <- matrix(parm[Data$pos.beta], Data$K, Data$T)
  zeta <- parm[Data$pos.zeta]
  parm[Data$pos.phi] <- phi <- interval(parm[Data$pos.phi], 1e-100, Inf)
  kappa <- interval(parm[Data$pos.kappa], 1e-100, Inf)
  parm[Data$pos.kappa] <- kappa
  lambda <- interval(parm[Data$pos.lambda], 1e-100, Inf)
  parm[Data$pos.lambda] <- lambda
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma
  parm[Data$pos.tau] <- tau <- interval(parm[Data$pos.tau], 1e-100, Inf)
  Sigma <- array(0, dim=c(Data$S, Data$S, Data$T))
  for (t in 1:Data$T) {
  Sigma[ , ,t] <- lambda[t]^2 * exp(-phi[t] * Data$D)^kappa[t]}
  ### Log-Prior
  beta.prior <- sum(dnormv(beta[,1], 0, 1000, log=TRUE),
  dnorm(beta[,-1], beta[,-Data$T], matrix(tau, Data$K,
  Data$T-1), log=TRUE))
  zeta.prior <- dmvn(zeta, rep(0,Data$S), Sigma[ , , 1], log=TRUE)
  phi.prior <- sum(dhalfnorm(phi[1], sqrt(1000), log=TRUE),
  dtrunc(phi[-1], "norm", a=0, b=Inf, mean=phi[-Data$T],
  sd=sigma[2], log=TRUE))
  kappa.prior <- sum(dhalfnorm(kappa[1], sqrt(1000), log=TRUE),
  dtrunc(kappa[-1], "norm", a=0, b=Inf, mean=kappa[-Data$T],
  sd=sigma[3], log=TRUE))
  lambda.prior <- sum(dhalfnorm(lambda[1], sqrt(1000), log=TRUE),
  dtrunc(lambda[-1], "norm", a=0, b=Inf, mean=lambda[-Data$T],
  sd=sigma[4], log=TRUE))
  sigma.prior <- sum(dhalfcauchy(sigma, 25, log=TRUE))
  tau.prior <- sum(dhalfcauchy(tau, 25, log=TRUE))
  ### Log-Likelihood
  mu <- tcrossprod(Data$X, t(beta))
  Theta <- matrix(zeta, Data$S, Data$T)
  for (t in 2:Data$T) {
  for (s in 1:Data$S) {
  Theta[s,t] <- sum(Sigma[,s,t] / sum(Sigma[,s,t]) * Theta[,t-1])}}
  mu <- mu + Theta
  LL <- sum(dnorm(Data$Y, mu, sigma[1], log=TRUE))
  ### Log-Posterior
  LP <- LL + beta.prior + zeta.prior + sum(phi.prior) +
  sum(kappa.prior) + sum(lambda.prior) + sigma.prior + tau.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
  yhat=rnorm(prod(dim(mu)), mu, sigma[1]), parm=parm)
  return(Modelout)
}

Initial.Values <- c(rep(0,S), rep(c(mean(Y),rep(0,K-1)),T), rep(1,T),
rep(1,T), rep(1,T), rep(1,4), rep(1,K))


Fit = LaplaceApproximation(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000 )
Initial.Values <- as.initial.values(Fit)
Fit = LaplacesDemon(Model, Data=Data, Initial.Values=Initial.Values, Iterations=1000, Status=100, Thinning=1 )
Fit = VariationalBayes(Model, Data=Data, parm=Initial.Values, Iterations=10000, Samples=1000, CPUs=5 )
Fit = IterativeQuadrature(Model, Data=Data, parm=Data$PGF(Data), Iterations=1000, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )

