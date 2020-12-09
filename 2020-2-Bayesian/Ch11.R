# Chapter11. JAGS 

library(rjags)


### 10.2 : mu0 = 10, sigma0^2 = 25 , a = 0.5 , b = 1 

# step2. Model Specification 
modelString = "
model{
  for(i in 1:n) {
  x[i] ~ dnorm(mu, invsigsq)
  }
  mu ~ dnorm(mu0, invsigsq0) 
  invsigsq ~ dgamma(a,b) 
  sigsq <- 1/invsigsq 
  mu0 <- 10 
  invsigsq0 <- 1/25 
  a <- 0.5 
  b <- 1 
}
"
writeLines(modelString, "model_ex12_1.txt")


# step3. data list 
dataList = list(n=10, x=c(10,13,15,11,9,18,20,17,23,21)) 

# step4. set initialization list 
initsList = list(mu=10, invsig=0.04)

# +) set init by resampling with replacement 
initsList = function(){
  resampledX = sample(x, replace=T) 
  muInit = sum(resampledX) / length(resampledX) 
  invsigsqInit = (1/sd(resampledX))^2 * 0.999 + 0.01 
  return(list(mu=muInit, invsigsq = invsigsqInit))
}

# step5. JAGS
jagsModel = jags.model(file="model_ex12_1.txt", data=dataList, inits=initsList, 
                       n.chains=3, n.adapt = 500)

# step6. burnin
update(jagsModel, n.iter=500)

# step7. MCMC samples 
codaSamples = coda.samples(jagsModel, variable.names = c("mu", "sigsq"), n.iter=5000)

# step8. traceplot 
coda::traceplot(codaSamples[,"mu"], main="", ylab="mu")
acf(codaSamples[,"mu"][[1]], plot=T, main="")

coda::traceplot(codaSamples[,"sigsq"], main="", ylab="sigsq")
acf(codaSamples[,"sigsq"][[1]], plot=T, main="")

# step9. posterior inference 
MuSamples = as.matrix(codaSamples[,"mu"])
SigSamples = as.matrix(codaSamples[, "sigsq"])

par(mfrow=c(1,2))
plot(density(MuSamples), xlab=expression(mu), ylab="posterior density", main="")
plot(density(SigSamples), xlab=expression(sigma^2), ylab="posterior density", main="")

# +) acceptance probability in JAGS 
AcceptRate = 1 - rejectionRate(codaSamples)
AcceptRate  # rejection rate = 0 
