# 2.3 Gibbs Sampling 

M <- 3000 ; m <- 500  # burn-in 500, nTotalIteration = 3000 
mu0 = 10 ; sigsq0 = 25 ; a=0.5 ; b=1 
x = c(10,13,15,11,9,18,20,17,23,21)

n = length(x) ; xbar <- mean(x) ; var.x <- var(x)  
THETA = matrix(nrow=M, ncol=2)  # empty list 
sigsq <- var.x ; mu <- xbar  # set initial value 

# Gibbs Sampler 
for (nsim in 1:M){
  
  # generate mu : Normal dist 
  condpost.var <- 1 / (1/sigsq0 + n/sigsq)
  condpost.mu  <- condpost.var * (xbar/(sigsq/n) + mu0/sigsq0)
  
  mu <- rnorm(1, condpost.mu, sqrt(condpost.var))
  
  # generate sigsq : Inverse Gamma dist 
  condpost.a <- a + n/2
  condpost.b <- b + 1/2*((n-1)*var.x + n*(xbar-mu)^2)
  
  sigsq <- 1/rgamma(1, condpost.a, condpost.b)
  
  # SAVE 
  THETA[nsim, ] <- c(mu, sigsq)
}

plot(x=THETA[,1] , y=THETA[,2], xlab="mu", ylab="sigma^2")    



