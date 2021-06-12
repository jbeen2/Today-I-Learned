# ======== 4.1 Gibbs Sampling ========

# 1. Linear Regression
# multivariate normal dist 
rmvnorm <- function(n, mu, Sig){  # n : #of data, mu, sig 
  p=length(mu) # p : dimension 
  R=chol(Sig)  # cholesky decomposition 
  z=matrix(rnorm(n*p),n,p)
  tt=z%*%R + matrix(mu, n, p, byrow=T)
  return(tt)
}

dat=read.csv("immigrants.csv")
y=dat$wage ; n=length(y)
X=cbind(rep(1,n),dat$sp, dat$lit) ; p=ncol(X)  # X design matrix 

a = 1 ; b = 1
XtX=t(X) %*% X
XtX.inv=solve(XtX) ; Xty=t(X)%*%y
beta_lse = as.vector(XtX.inv %*% t(X) %*% y)  # beta hat 
sigsq.hat=sum((y - X %*% beta_lse)^2)/(n - p) # sigsq hat : 1.0866 

# initial value setting 
beta0=beta_lse
Sig0=diag(diag(XtX.inv))*sigsq.hat*100  
# > diag(diag(XtX.inv)) = diag matrix
# > sigsq * XtX : correlation 있는 matrix, diag element 만 가져다가 쓸 것 
# > * 100 = prior 정보 넉넉히 사용할 것 
Sig0.inv=solve(Sig0)

N = 10000 ; nburn=1000  # #of iterations 
sigsq.samples = rep(0,N) ; beta.samples = matrix(0, N, p) # 값을 저장할 empty matrix 생성 
# start point 
beta.init=beta_lse
sigsq.init=sigsq.hat
beta=beta.init ; sigsq=sigsq.init

# start gibbs sampling
for (iter in 1:(N+nburn)) {
  Sig.beta = solve(Sig0.inv + XtX/sigsq)
  mu.beta = Sig.beta %*%(Sig0.inv%*%beta0 + 1/sigsq*Xty)
  beta = as.vector(rmvnorm(1, mu.beta, Sig.beta))
  SSR = sum((y - X %*% beta)^2)/(n - p)
  sigsq = 1/rgamma(1, n/2 + a, 1/2 *SSR + b)
  if (iter > nburn) {
    beta.samples[iter-nburn,] = beta
    sigsq.samples[iter-nburn] = sigsq
}
}

# 95% HPD interval : Trace Plot 
ci_beta = round(apply(beta.samples, 2, quantile, probs = c(0.025, 0.975)),4)
ci_sigsq = round(quantile(sigsq.samples, c(0.025, 0.975)),4)
par(mfrow=c(2,2))
for (i in 1:3) plot(beta.samples[,i], type="l",xlab=paste0("beta_",i), ylab="", col="blue")
plot(sigsq.samples,xlab="sigsq",ylab="", type="l", col="blue")
# > Thinning 을 하지 않아서 뭉치는 경향이 있지만, good! 

# 모수의 사후분포 & 95% 사후구간 
par(mfrow=c(2,2))
for (i in 1:3) {
  plot(density(beta.samples[,i]),
       type="l",main="",xlab=paste0("beta_",i))
  abline(v=ci_beta[,i], col=2, lty=2) }
plot(density(sigsq.samples),main="",xlab="sigsq", type="l")
abline(v=ci_sigsq, col=2, lty=2)


# 2. Truncated Multivate Normal Distribution 
library(truncnorm) # univariate normal 
k = 5 ; mu = c(0,1,2,3,5) ; Sig = matrix(0.7, 5, 5) + diag(k) * 0.3 ; invSig = solve(Sig)
m = 1000 ; N = 10000  # m : burnin, M : iteration 
theta.init = c(0,1,2,3,4)  # 제한조건을 만족하는 초기치 선택 
theta.samples = matrix(0, N, k)
theta=c(1:k)

for (iter in 1:(m+N)) {
  for (i in 1:k) {
    vec.Ai = invSig[i,-i]
    vec.mi = (theta-mu)[-i]
    cond.mean = mu[i] - 1/invSig[i,i] * vec.Ai %*% vec.mi
    cond.sd = 1/sqrt(invSig[i,i])
    a = ifelse(i==1, Inf, theta[i-1]) # -inf < theta1 < theta2 
    b = ifelse(i==k, Inf, theta[i+1]) # upper bound inf 
    theta[i] = rtruncnorm(1, a,b, cond.mean, cond.sd)
  }
  if (iter > m) theta.samples[iter-m,] = theta
}

# 모수의 산점도 
