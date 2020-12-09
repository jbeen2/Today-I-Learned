library(rjags)

modelString = "
model{
for (i in 1:length(y)){ 
y[i] ~ dnorm(inprod(x[i,], beta[]), invsigsq)
}
beta[1:length(beta0)] ~ dmnorm(beta0[], invSig0[,])

invsigsq ~ dgamma(a,b)
sigsq = 1/invsigsq

invSig0 = 1/(c*sigsq) * xtx
}
"
writeLines(modelString, "model_final.txt")

# read data 
icecream = read.table("icecream.txt", header=T)
y = icecream[,2] 
x1 = icecream[,3] ; x2 = icecream[,4] ; x3 = icecream[,5]
n = length(y)
p=4
K=p
x = matrix(0,n,p)
# x = icecream[,3:5]
for (i in 1:n) {
  x[i,1] = 1
  x[i,2] = x1[i] 
  x[i,3] = x2[i]
  x[i,4] = x3[i]
}
xtx = t(x)%*%x
xtx.inv = solve(xtx)
xty = t(x) %*% y
beta.LSE = as.vector(xtx.inv %*% xty)
SSR = t(y) %*% (diag(1,nrow=n) - x %*% xtx.inv %*% t(x) ) %*% y
sigsq.hat = as.numeric(SSR) / n-p

# prior 
beta0 = beta.LSE
c=n ; a = 0.5 ; b = a*sigsq.hat

dataList = list(x=x, y=y, xtx=xtx, a=a, b=b, c=c, beta0=beta0)
initsList = list(beta=beta0, invsigsq = 1/sigsq.hat)

nChains=3

jagsModel = jags.model(file="model_final.txt", data=dataList, 
                       inits = initsList, n.chains = nChains, n.adapt=100)

update(jagsModel, n.iter=1000)
codaSamples = coda.samples(jagsModel, variable.names = c("beta", "sigsq"), n.iter=30000)
variable.names(codaSamples[[1]])

para = c("beta0", "beta1", "beta2", "beta3")
par(mfrow=c(2,4))
for(i in 1:4){
  coda::traceplot(codaSamples[,i], main="", ylab=para[i])
  acf(codaSamples[,i][[1]], plot=T, main="")
}

MCMCsamples = as.matrix(codaSamples)
par(mfrow=c(2,2))
for(i in 1:4) plot(density(MCMCsamples[,i]), main="", xlab=para[i])

# beta0 
beta0.sample = as.matrix(codaSamples[, "beta[1]"])
mean(beta0.sample) ; sd(beta0.sample)  # 153.3477, 141.627 

# beta1 
beta1.sample = as.matrix(codaSamples[, "beta[2]"])
mean(beta1.sample) ; sd(beta1.sample)  # -0.53, 0.403 

# beta2
beta2.sample = as.matrix(codaSamples[, "beta[3]"])
mean(beta2.sample) ; sd(beta2.sample)  # 0.00042, 0.00014

# beta3
beta3.sample = as.matrix(codaSamples[, "beta[4]"])
mean(beta3.sample) ; sd(beta3.sample)  # 3.4898, 0.4244
