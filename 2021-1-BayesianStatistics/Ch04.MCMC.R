mu0 <- 10 ; sig0sq <- 25 ; a <- .5 ; b <- 1 
x <- c(10,13,15,11,9,18,20,17,23,21)
dataList = list(x=x, mu0=mu0, sig0sq=sig0sq, a=a, b=b)

# posterior kernel 
post.normal_mu_sigsq = function(theta, dataList){
  x = dataList$x ; mu0 = dataList$mu0 ; sig0sq = dataList$sig0sq ; 
  a = dataList$a ; b = dataList$b 
  
  mu = theta[1] ; sigsq = theta[2] 
  
  f = exp(- 0.5 * length(x) * log(sigsq) - 0.5 * sum((x-mu)^2) / sigsq 
          - 0.5 * (mu-mu0)^2 / sigsq
          - (a+1) * log(sigsq) - b/sigsq)
  
  return (f)
}

Metropolis_normal_mu_sigsq = function(nsim, nburn, delta, dataList, initsList){ 
  mu = initsList$mu ; log.sigsq = log(initsList$sigsq) 
  theta.curr = c(mu, log.sigsq)
  p = length(theta.curr)
  
  para.samples = matrix(0, nsim, p) 
  for (iter in 1:(nsim+nburn)){ 
    z = rnorm(p, 0, 1) 
    theta.prop = z * delta + theta.curr 
    mu.curr = theta.curr[1] 
    sigsq.curr = exp(theta.curr[2]) 
    mu.prop = theta.prop[1] 
    sigsq.prop = exp(theta.prop[2]) 
    alpha = post.normal_mu_sigsq(c(mu.prop, sigsq.prop), dataList) / post.normal_mu_sigsq(c(mu.curr, sigsq.curr), dataList) * sigsq.prop / sigsq.curr
    
    if (runif(1) < alpha) {theta.next <- theta.prop} else {theta.next <- theta.curr}
    theta.curr = theta.next 
    
    if (iter > nburn) para.samples[iter-nburn, ] = c(theta.next[1] , exp(theta.next[2]))
    }
  return(para.samples)
}

nChains = 3 
nsim = 20000 ; nburn=5000 ; p=2 
mcmc.samples = array(0, dim=c(nsim, p, nChains))

delta = 1 

inits.random = function(x){
  resampledX = sample(x, replace=T) 
  muInit = mean(resampledX) ; sigsqInit = var(resampledX)
  return(list(mu = muInit, sigsq = sigsqInit))
}

for (ich in 1:nChains) { 
  initsList = inits.random(x)
  mcmc.samples[,,ich] = Metropolis_normal_mu_sigsq(nsim, nburn, delta, dataList, initsList)
}


### posterior inference 
mcmc.samples.combined = rbind(mcmc.samples[,,1], mcmc.samples[,,2], mcmc.samples[,,3])
para.hat = apply(mcmc.samples.combined, 2, mean)
HPD = apply(mcmc.samples.combined, 2, function(x) quantile(x, c(0.025, 0.975)))

para.hat


##### (2) Multivariate MCMC 
post.normal_mu_sigsq2 = function(theta, dataList){
  x = dataList$x ; mu0 = dataList$mu0 ; sig0sq = dataList$sig0sq ; 
  a = dataList$a ; b = dataList$b 
  
  mu = theta[1] ; sigsq1 = theta[2] ; sigsq2 = theta[3]
  sigsq = sqrt(sigsq1 * sigsq2)
  
  f = exp(- 0.5 * length(x) * log(sigsq) - 0.5 * sum((x-mu)^2) / sigsq 
          - 0.5 * (mu-mu0)^2 / sigsq
          - (a+1) * log(sigsq) - b/sigsq)
  
  return (f)
}


library(truncnorm)
Metropolis_normal_mu_sigsq2 = function(nsim, nburn, delta, dataList, initsList){ 
  mu = initsList$mu ; log.sigsq1 = log(initsList$sigsq) ; log.sigsq2 = log(initsList$sigsq) 
  theta.curr = c(mu, log.sigsq1, log.sigsq2)
  p = length(theta.curr)
  para.samples = matrix(0, nsim, p) 

  for (iter in 1:(nsim+nburn)){ 
    z = rtruncnorm(p, 0, 1) 
    theta.prop = z * delta + theta.curr 
    mu.curr = theta.curr[1] 
    sigsq1.curr = exp(theta.curr[2]) 
    sigsq2.curr = exp(theta.curr[3]) 
    
    mu.prop = theta.prop[1] 
    sigsq1.prop = exp(theta.prop[2]) 
    sigsq2.prop = exp(theta.prop[3]) 
    
    alpha = post.normal_mu_sigsq2(c(mu.prop, sigsq1.prop, sigsq2.prop), dataList) / post.normal_mu_sigsq2(c(mu.curr, sigsq1.curr, sigsq2.curr), dataList) 
    
    if (runif(1) < alpha) {theta.next <- theta.prop} else {theta.next <- theta.curr}
    theta.curr = theta.next 
    
    if (iter > nburn) para.samples[iter-nburn, ] = c(theta.next[1] , exp(theta.next[2]), exp(theta.next[3]))
  }
  return(para.samples)
}

nChains = 3 
nsim = 20000 ; nburn=5000 ; p=3
mcmc.samples = array(0, dim=c(nsim, p, nChains))

delta = 1 

inits.random = function(x){
  resampledX = sample(x, replace=T) 
  muInit = mean(resampledX) ; sigsqInit = var(resampledX)
  return(list(mu = muInit, sigsq = sigsqInit))
}

for (ich in 1:nChains) { 
  initsList = list(mu = para.hat[1], sigsq = para.hat[2])
  mcmc.samples[,,ich] = Metropolis_normal_mu_sigsq2(nsim, nburn, delta, dataList, initsList)
}


# (4) Acceptance Rate 
require(coda)
Metro.draws = mcmc(mcmc.samples[,,1])
accept.rate = 1 - rejectionRate(Metro.draws)
accept.rate


mu.samples = mcmc.samples[,1,] ; sigsq.samples = mcmc.samples[,2,]
par(mfrow=c(2,2))
plot(mu.samples[,1], type="l", xlab="iteration", ylab=quote(mu),
     main = paste0("accept.rate = ", round(accept.rate[1],3)))
lines(density(mu.samples[,2]), col=2)
lines(density(mu.samples[,3]), col=3)
