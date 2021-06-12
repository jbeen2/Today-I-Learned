setwd("~/workdir/statistics/Today-I-Learned/2021-1-BayesianStatistics")

library(MASS)
library(dplyr)
library(tidyverse)
library(rjags)
library(data.table)

# Used Car Price Prediction 
data <- read.csv("UsedCar2.csv")

min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
data_norm <- as.data.frame(lapply(data %>% select(-log_Price), min_max_norm))
data_norm$Price <- data$log_Price


# ======================================
# 1. linear regression 
# ======================================
result <- lm(Price ~ . , data=data_norm)
summary(result)


# Stepwise 
stepAIC(result, direction = 'both')



# ======================================
# 2. GVS (Gibbs Variable Selection) 
# pseudo-prior : Normal dist, mean/var = LSE 
# ======================================
mu.beta <- result$coefficients
var.beta <- diag(vcov(result))

y <- data_norm$Price ; x <- data_norm %>% select(-Price) ; n <- nrow(data)
X <- cbind(rep(1,n), x) ; K <- ncol(X)

# pseudo prior 
pseudo.beta.mean <- mu.beta
pseudo.beta.var <- var.beta

# prior 
prior.beta.var <- var.beta * 100
prior.beta.mean <- rep(0, K)


# Model Setting for JAGS 
modelString = "
model {
  for(j in 1:K){gbeta[j] <- gamma[j]*beta[j]}
  for(i in 1:n){
    y[i] ~ dnorm(mu[i], invsigsq)
    mu[i] <- inprod(X[i,1:K], gbeta [1:K])
  }
  
  for(j in 1:K){gamma[j] ~ dbern (0.5)}
  for(j in 1:K){
      beta[j] ~ dnorm (m.b[j], tau.b[j])
      m.b[j] <- gamma[j] * prior.beta.mean[j] + (1+gamma[j]) * pseudo.beta.mean[j]
      tau.b[j] <- gamma[j] / prior.beta.var[j] + (1+gamma[j]) / pseudo.beta.var[j]
  }
    
  invsigsq ~ dgamma (0.01, 0.01)
  }
  "
write(modelString , file="model_GVS2.txt")


# MCMC 
dataList <- list(n=n, K=K, y=y, X=X, 
                 pseudo.beta.mean=pseudo.beta.mean, pseudo.beta.var=pseudo.beta.var, 
                 prior.beta.var=prior.beta.var, prior.beta.mean=prior.beta.mean)
gammaInit=rep(1,K) # gamma init = 1 
initsList=list(beta=mu.beta, gamma=gammaInit) 
            
nChains=3; nIter=10000
jagsModel=jags.model(file="model_GVS2.txt", data=dataList, inits=initsList,
                     n.chains=nChains, n.adapt = 1000)
update(jagsModel, n.iter = 3000)
codaSamples=coda.samples(jagsModel, variable.names=c("gamma", "beta"),
                         n.chains=nChains, n.iter=nIter)
para.samples=as.matrix(codaSamples)

head(para.samples)
beta.samples = para.samples[, 1:K] # 1 ~ 14 
gamma.samples = para.samples[, (K+1) : (K+K)] # 15 ~ 28 


# Variable Selection 
m = gamma.samples
mm = as.data.table(m)[, .N, by=eval(paste0("gamma[", seq_len(ncol(m)), "]"))]
mm.order = order(mm$N, decreasing=T)
mm$N = round(mm$N/(nIter*nChains), 4)
gamma.hat = as.numeric(mm[which.max(mm$N)])
gamma.hat = gamma.hat[1:K] ; gamma.hat

mm[mm.order[1:10]]


# Model1 Coef Estimation 
gamma.samples.collapsed <- apply(gamma.samples, 1, function(x) paste(x, collapse=" "))
gamma.hat.collapsed = paste(gamma.hat , collapse=" ")
id.selected = which(gamma.samples.collapsed == gamma.hat.collapsed)
length(id.selected) # 17750 = 3 * 10000 * 0.5917 

beta.samples.selected = beta.samples[id.selected, ]
colnames(beta.samples.selected)=c("intercept", colnames(x))
beta.samples.selected2 = beta.samples.selected[, gamma.hat == 1]
head(beta.samples.selected2)
beta.selected.hat = apply(beta.samples.selected2, 2, function(x) quantile(x, c(0.025,0.5,0.975)))
t(beta.selected.hat)



# ======================================
# 3. Bayesian Linear Regression 
# ======================================
y <- data_norm$Price ; x2 <- data_norm[,colnames(beta.selected.hat)] ; n <- nrow(data_norm)
# X2 <- cbind(rep(1,n), x2) ; # without intercept 
K2 <- ncol(x2)

a = 1 ; b = 1 ; X2 <- as.matrix(x2)
XtX = t(X2)%*%X2 ; XtX.inv = solve(XtX) ; Xty = t(X2)%*%y
beta.hat <- as.vector(XtX.inv %*% Xty)
sigsq.hat <- sum((y-X2 %*% beta.hat)^2)/(n-K2)
beta0 <- beta.hat 
Sig0 <- diag(diag(XtX.inv))*sigsq.hat*100
Sig0.inv <- solve(Sig0)

modelString2 = "
model
{
  for (i in 1:length(y)){
    y[i] ~ dnorm(inprod(X[i,], beta[]), invsigsq)
  }
  beta[1:length(beta0)] ~ dmnorm(beta0[], Sig0.inv[,])
  
  invsigsq ~ dgamma(a,b)
  sigsq = 1/invsigsq
}
"
writeLines(modelString2, "model_reg.txt")

dataList2 = list(X=X2, y=y, a=a, b=b, beta0=beta0, Sig0.inv=Sig0.inv)
initsList2 = list(beta=beta.hat, invsigsq = 1/sigsq.hat)
nChains=3

jagsModel2 = jags.model(file="model_reg.txt", data=dataList2, inits=initsList2, 
                        n.chains = nChains, n.adapt = 3000)
update(jagsModel2, n.iter=1000)
codaSamples2=coda.samples(jagsModel2, variable.names=c("beta", "sigsq"), n.iter=10000) # 30000
para.names2 = variable.names(codaSamples2[[1]])


# Trace Plot & ACF 
par("mar")
par(mar = c(1,1,1,1))
par(mfrow=c(4,2))
for(i in 1:8) { 
  coda::traceplot(codaSamples2[,i], main="", ylab=para.names2[i])
  acf(codaSamples2[,i][[1]], plot=T, main=para.names2[i])
}

MCMCSamples = as.matrix(codaSamples2)
HPD = round(apply(MCMCSamples, 2, quantile, probs=c(.025, .975)), 4)


# Posterior density & 95% HPD
par(mfrow=c(4,2))
for(i in 1:8) {
  plot(density(MCMCSamples[,i]), main="", xlab=para.names2[i], col='blue')
  abline(v=HPD[,i], col=2)
}


# Gelman Statistics 
gelman = gelman.diag(codaSamples2) # beta 
gelman.1 = as.matrix(gelman$psrf)
if (max(gelman.1 > 1.1)) cat("Warning : Gelman Shrink Factor > 1.1", "\n")
gelman.2 = gelman$psrf
if (max(gelman.2 > 1.1)) cat("Warning : Gelman Multivariate Shrink Factor > 1.1", "\n")
gelman.1


# Acceptance Rate 
AcceptRate = 1 - rejectionRate(codaSamples2) ; AcceptRate


# beta hat 
mcmc.beta.hat = apply(MCMCSamples, 2, function(x) quantile(x, c(0.025,0.5,0.975)))
colnames(mcmc.beta.hat)=c(colnames(X2), "sigsq")
t(mcmc.beta.hat)


