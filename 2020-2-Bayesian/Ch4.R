# Ch.4 Bayesian Inference 


# Highest Posterior Density (HPD) interval 

# 2) using grid point 
HPDgrid <- function(prob, level=0.95){
  prob.sort = sort(prob, decreasing=T)
  M = min(which(cumsum(prob.sort) >= level))
  height = prob.sort[M]
  HPD.index = which(prob>=height)
  HPD.level = sum(prob[HPD.index])
  res = list(index = HPD.index, level=HPD.level)
  return(res)
}


N=1001 ; theta = seq(-3,3,length=N)
prob = exp(-.5 / .25*(theta-.3)**2)
prob = prob / sum(prob)
alpha = 0.05 ; level = 1-alpha


HPD = HPDgrid(prob, level) ; HPD
plot(theta[HPD$index] , rep(1, length(HPD$index)))
HPDgrid.hat <- c(min(theta[HPD$index]), max(theta[HPD$index]))
HPDgrid.hat ; HPD$level



# 3) using posterior samples 
HPDsample = function(theta, level=0.95) { 
  N = length(theta) ; theta.sort = sort(theta)
  M = ceiling(N*level) 
  nCI = N-M   # number of possible CIs 
  CI.width = rep(0, nCI)
  
  for (i in 1:nCI) CI.width[i] = theta.sort[i+M] - theta.sort[i]
  index = which.min(CI.width)
  HPD = c(theta.sort[index], theta.sort[index+M])
  return(HPD)
  }

N = 30000 ; theta = rnorm(N, 0.3, 0.5)
HPDsample(theta)




# 4 
N=1001 ; x = seq(0, 10, length=N)
prob = dbeta(x, 4, 8)
prob = prob / sum(prob)

alpha = 0.05 ; level = 1-alpha

HPD = HPDgrid(prob, level) ; HPD
plot(x[HPD$index] , rep(1, length(HPD$index)))
HPDgrid.hat <- c(min(x[HPD$index]), max(x[HPD$index]))
HPDgrid.hat ; HPD$level


# 5 
N = 100 ; theta = rnorm(N, 14.97, 3.93)
HPDsample(theta)


# 6 
theta <- seq(0.1, 0.9, .1) ; x=5 ; n=20 
joint <- theta**x * (1-theta)**(n-x)
post <- joint / sum(joint)
post

sum(post[1:4])
