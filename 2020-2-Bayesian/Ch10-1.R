# MCMC 

### 10.2 Gibbs Sampling 

# initial values 
M <- 3000 ; m <- 500 # burnin size 
mu0 <- 10 ; sig20 <- 25 # mu ~ N(10, 5^2)
a <- 0.5 ; b<- 1 # sigma^2 ~ IG(0.5, 1)
x <- c(10,13,15,11,9,18,20,17,23,21)
n <- length(x) 
xbar <- mean(x) ; var.x <- var(x) 
THETA <- matrix(nrow=M , ncol=2) # theta sample 생성 

# initial value of sig2 
sig2 <- var.x # 표본분산으로 초기치 잡는다 
condpost.a <- a+n/2 

# simulation 
for (nsim in 1:M) {
  
  # generate mu : post(pi) ~ N (7장)
  condpost.mu <- (sig2/n * mu0 + sig20 * xbar) / (sig2/n + sig20) 
  condpost.var <- 1/(1/sig20 + n/sig2) 
  mu <- rnorm(1, condpost.mu, sqrt(condpost.var))  # 초기치 하나 생성됨 
  
  # generate sig2 
  # 위에서 만들어 낸 초기치 사용 
  condpost.b <- b+1/2 * ((n-1)*var.x + n*(xbar-mu)^2) 
  sig2 <- 1/rgamma(1, condpost.a, condpost.b)  # IG 
  
  # save 
  THETA[nsim, ] <- c(mu, sig2)
  
}

THETA



### 10.3 Gibbs Sampling 
# y ~ Multi, prior 
# Y : observed <-> X : not observed 

Y <- c(14,1,1,1,5) ; n<-22 ; 
a<-c(.25, .25, .25, .25) ; b <- c(.125, 0,0,0.375)  # a1 ~ a4, b1 ~ b4 
c<- 0.5 # Y5 
alpha <- c(1,1,1) ; m<- 1000 ; M <- 3000 
THETA <- matrix(0, m+M , 2) ; X <- matrix(0, m+M, 9)

# Dirichelet generator 
rDirich <- function(n,p,a){
  rDir <- matrix(0,n,p)
  for (nsim in 1:n) { 
    x <- c(1: (p+1))*0  # np.zeros(1,p+1)
    for (i in 1:(p+1)) x[i] <- rgamma(1,a[i])
    s <- sum(x) 
    for ( i in 1:(p+1)) x[i] <- x[i] / s 
    rDir[nsim, 1:p] <- x[1:p]
  }
  rDir
}

# initialize 
theta <- 0.5 ; eta <- 0.5 
THETA[1, 1:2] <- c(theta, eta) 
X[1,1] <- 0.5*Y[1] ; X[1,3] <- 0.5*Y[2] ; X[1,5] <- 0.5*Y[3] 
X[1,7] <- 0.5*Y[4] ; X[1,9] <- Y[5]

# simulate 
for (i in 2:(m+M)) { 
  p <- 2 
  aa <- c(X[i-1,1] + X[i-1,3] + alpha[1], 
          X[i-1,5] + X[i-1,7] + alpha[2], 
          X[i-1,9] + alpha[3])
  
  THETA[i,] <- rDirich(1,p,aa)
  theta <- THETA[i,1] ; eta<-THETA[i,2] 
  
  X[i,1] <- rbinom(1,Y[1], a[1]*theta / (a[1]*theta+b[1]))
  X[i,3] <- rbinom(1,Y[2], a[2]*theta / (a[2]*theta+b[2]))
  X[i,5] <- rbinom(1,Y[3], a[3]*eta / (a[3]*eta+b[3]))
  X[i,7] <- rbinom(1,Y[4], a[4]*eta / (a[4]*eta+b[4]))
  
  X[i,2] <- Y[1] - X[i,1] 
  X[i,4] <- Y[2] - X[i,3] ; X[i,6] <- Y[3] - X[i,5]
  X[i,8] <- Y[4] - X[i,7] ; X[i,9] <- Y[5]
  
}

theta.sim <- THETA[(m+1):(m+M) , 1]
eta.sim <- THETA[(m+1):(m+M) , 2]

# end simulation 
par(mfrow=c(1,1)) 
plot(theta.sim, eta.sim, xlab=expression(theta), ylab=expression(eta))

par(mfrow=c(1,2))
plot(density(theta.sim), type="l", xlab=expression(theta), ylab='posterior', main='')
plot(density(eta.sim), type="l", xlab=expression(eta), main='')