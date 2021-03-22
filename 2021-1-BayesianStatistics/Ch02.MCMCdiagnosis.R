# MCMC convergence & diagnosis 

library(coda)
n = 20 ; xbar = 4 ; sigma = 1 
mu0 = 0 ; sigma0 = 10  
dataList = list(n=n, xbar=xbar, sigma=sigma, mu0=mu0, sigma0=sigma0)

# Function to Compute Posterior Kernel 
post.kernel = function(mu, dataList){

  post.kernel = exp(-.5 * (((dataList$xbar - mu) * sqrt(dataList$n) / dataList$sigma)^2 
                           + ((mu - dataList$mu0) / dataList$sigma0)^2 ))
  return(post.kernel)
}

# Perform Random Walk Metropolis 
Metro = function(nsim, mu.init, delta, dataList){
  mu.samples = mu.init 
  mu.curr = mu.init 
  
  for (iter in 1:nsim){
    mu.prop = rnorm(1, mean=mu.curr, sd=delta)
    
    # alpha = (proposal | x) / (current | x) 
    alpha = min(1, post.kernel(mu.prop, dataList) / post.kernel(mu.curr, dataList))
    
    # generate random value  
    u = runif(1) 
    
    # u<alpha -> mu.proposal & u>alpha -> mu.current 
    mu.next = mu.prop * (u<alpha) + mu.curr * (u>alpha)
    mu.samples = rbind(mu.samples, mu.next)
    
    # 현재 반복에서의 next = 다음 반복에서의 current 
    mu.curr = mu.next
  }
  return (mu.samples)
}

######################################## 

delta = .2 ; nsim = 10000 ; n.chains = 3 
mu.Samples = matrix(0, nsim, n.chains)
for (i in 1:n.chains) {
  mu.init = rnorm(1, mu0, 2)
  mu.Samples[, i] = Metro(nsim-1, mu.init, delta, dataList)
}

# trace plots : 501 ~ 10000 
nwarm = 500 ; nsim = 10000 
plot(mu.Samples[(nwarm+1) : nsim, 1], type='l')
lines(mu.Samples[(nwarm+1) : nsim, 2], col=2)
lines(mu.Samples[(nwarm+1) : nsim, 3], col=3)

# posterior pdf plots of mu : 1 ~ 200 
nwarm = 0 ; nsim = 200
plot(mu.Samples[(nwarm+1) : nsim, 1], type='l')
lines(mu.Samples[(nwarm+1) : nsim, 2], col=2)
lines(mu.Samples[(nwarm+1) : nsim, 3], col=3)

# plot of Gelman psf (potential scale reduction factor)
mu3.Samples = mcmc(mu.Samples[2001:10000, ]) # samples of mu from 3 chains 
mu3.codaSamples = mcmc.list(list(mu3.Samples[,1], mu3.Samples[,2], mu3.Samples[,3]))
gelman.plot(mu3.codaSamples, col=c("black", "blue"))

# ACF plot of mu[1] for chain 1 
plot(acf(mu.Samples[(nwarm+1):nsim, 1]), main = "Autocorrelation")

# ACFs 
mu.postSamples = as.matrix(mu3.codaSamples)
aa = acf(mu.postSamples) ; aa
