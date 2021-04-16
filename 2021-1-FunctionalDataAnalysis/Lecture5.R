##### Lecture 5 Part 1 #####

library(fda)
library(refund)
set.seed(2016)

par(mfrow=c(2,2),mar=c(2,2,1,1))
gaus2<-function(x){exp(-x^2/2)/sqrt(2*pi)}
gaus4<-function(x){(1/2)*(3-x^2)*gaus2(x)}
ep2<-function(x){(3/4)*(1-x^2)*(abs(x)<1)}
ep4<-function(x){(15/8)*(1-(7/3)*x^2)*ep2(x)}

xlim=4
ylimh=1.4
yliml=-0.25
fnt = 1.3

par(mfrow=c(2,2),mar=c(2,4,1,1), cex=fnt)
curve(gaus2,from=-xlim,to=xlim,xlab="",ylab="Gaussian - 2nd",ylim=c(yliml,ylimh))
abline(h = 0,lty=2)
curve(gaus4,from=-xlim,to=xlim,xlab="",ylab="Gaussian - 4th",ylim=c(yliml,ylimh))
abline(h = 0,lty=2)
curve(ep2,from=-xlim,to=xlim,xlab="",ylab="Epanechnikov - 2nd",ylim=c(yliml,ylimh))
abline(h = 0,lty=2)
curve(ep4,from=-xlim,to=xlim,xlab="",ylab="Epanechnikov - 4th",ylim=c(yliml,ylimh))
abline(h = 0,lty=2)

N = 100; M = 5
T = matrix(runif(N*M),nrow=N,ncol=M)
T = apply(T,1,"sort"); T = t(T)

mu <-function(t){
  .5*dnorm(t,mean = 1/3,sd=1/18)+
    .1*dnorm(t,mean = 2/3,sd=1/18)}


par(mar=c(2,2,1,1))
curve(mu,from=0,to=1)

library(MASS)
C<-function(d){exp(-abs(d))}
Y<-matrix(0,nrow=N,ncol=M)
for(n in 1:N){
  tms<-T[n,]
  Sig<-C(outer(tms,tms,FUN="-"))
  mu_v<-mu(tms)
  Y[n,] = mvrnorm(1,mu_v,Sig)
}

par(mar=c(2,2,1,1))
plot(c(T),c(Y),xlab="",ylab="")


h = 1/10; P<-1; grd_n<-50
K<-function(d){(1/sqrt(2*pi))*exp(-d^2/2)}
grd<-seq(0,1,length=grd_n)
mu_hat = numeric(grd_n)
y = c(Y); t_obs = c(T)
for(i in 1:grd_n){
  t = grd[i]
  x = matrix(0,nrow=N*M,ncol=P+1)
  for(j in 0:P){
    x[,j+1] = (t-t_obs)^j 
  }
  W = diag(K((t-t_obs)/h))
  beta_hat = solve(t(x)%*%W%*%x,t(x)%*%W%*%y)
  mu_hat[i] = beta_hat[1]
}

par(mar=c(2,2,1,1))
plot(c(T),c(Y),xlab="",ylab="")
points(grd,mu_hat,type="l",lwd=1.5,col="red")

newdata=data.frame(t_obs = grd)
fit1<-loess(y~t_obs,span=1/2)
pred1 <- predict(fit1,newdata=newdata,
                 degree=1,family="gaussian")
fit2<-loess(y~t_obs,span=1/10)
pred2 <- predict(fit2,newdata=newdata,
                 degree=1,family="gaussian")
fit3<-loess(y~t_obs,span=1/15)
pred3 <- predict(fit3,newdata=newdata,
                 degree=1,family="gaussian")
fit4<-loess(y~t_obs,span=1/20)
pred4 <- predict(fit4,newdata=newdata,
                 degree=1,family="gaussian")


par(mar=c(2,2,1,1))
plot(t_obs,y,xlab="",ylab="")
points(grd,pred1,type="l",col=2,lwd=2)
points(grd,pred2,type="l",col=3,lwd=2)
points(grd,pred3,type="l",col=4,lwd=2)
points(grd,pred4,type="l",col=5,lwd=2)


##### Lecture 5 Part 2 #####

library(fda)
library(refund)
set.seed(2016)


  N = 100; M = 5
T = matrix(runif(N*M),nrow=N,ncol=M)
T = apply(T,1,"sort"); T = t(T)
mu <-function(t){
  .5*dnorm(t,mean = 1/3,sd=1/18)+
    .1*dnorm(t,mean = 2/3,sd=1/18)}
library(MASS)
C<-function(d){exp(-abs(d))}
Y<-matrix(0,nrow=N,ncol=M)
for(n in 1:N){
  tms<-T[n,]
  Sig<-C(outer(tms,tms,FUN="-"))
  mu_v<-mu(tms)
  Y[n,] = mvrnorm(1,mu_v,Sig)}

par(mar=c(2,2,1,1))
plot(c(T),c(Y),xlab="",ylab="")

y<-c(Y); t<-c(T)
K<-exp(-outer(t,t,FUN="-")^2)
lambda1<-1
alpha1<-solve(lambda1*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})

grd<-seq(0,1,length=50)
K_fit<-exp(-outer(grd,t,FUN="-")^2)
muhat<-K_fit%*%alpha1

par(mar=c(2,2,1,1))
plot(c(T),c(Y),xlab="",ylab="")
points(grd,muhat,type="l",lwd=2)


y<-c(Y); t<-c(T)
K<-exp(-outer(t,t,FUN="-")^2)
lambda1<-1
alpha1<-solve(lambda1*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda2<-.01
alpha2<-solve(lambda2*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda3<-.0001
alpha3<-solve(lambda3*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda4<-.000001
alpha4<-solve(lambda4*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})


grd<-seq(0,1,length=50)
K_fit<-exp(-outer(grd,t,FUN="-")^2)
muhat1<-K_fit%*%alpha1
muhat2<-K_fit%*%alpha2
muhat3<-K_fit%*%alpha3
muhat4<-K_fit%*%alpha4


par(mar=c(2,2,1,1))
plot(c(T),c(Y),xlab="",ylab="")
points(grd,muhat1,type="l",lwd=2,col=1)
points(grd,muhat2,type="l",lwd=2,col=2)
points(grd,muhat3,type="l",lwd=2,col=3)
points(grd,muhat4,type="l",lwd=2,col=4)
leg<-paste("lambda=",c(1,.01,.0001,.000001))
legend("topright",legend=leg,col=1:4)


y<-c(Y); t<-c(T)
K<-exp(-abs(outer(t,t,FUN="-")))
lambda1<-1
alpha1<-solve(lambda1*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda2<-.5
alpha2<-solve(lambda2*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda3<-.1
alpha3<-solve(lambda3*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda4<-.05
alpha4<-solve(lambda4*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})


grd<-seq(0,1,length=100)
K_fit<-exp(-abs(outer(grd,t,FUN="-")))
muhat1<-K_fit%*%alpha1
muhat2<-K_fit%*%alpha2
muhat3<-K_fit%*%alpha3
muhat4<-K_fit%*%alpha4


par(mar=c(2,2,1,1))
plot(c(T),c(Y),xlab="",ylab="")
points(grd,muhat1,type="l",lwd=2,col=1)
points(grd,muhat2,type="l",lwd=2,col=2)
points(grd,muhat3,type="l",lwd=2,col=3)
points(grd,muhat4,type="l",lwd=2,col=4)
leg<-paste("lambda=",c(1,.5,.1,.05))
legend("topright",legend=leg,col=1:4)


y<-c(Y); t<-c(T)
myfun<-function(t,s){
  if(t<=s){cosh(1-s)*cosh(t)
  }else{cosh(1-t)*cosh(s)}
}
myfun<-Vectorize(myfun)
K<-outer(t,t,FUN=myfun)

  

lambda1<-1
alpha1<-solve(lambda1*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda2<-.5
alpha2<-solve(lambda2*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda3<-.1
alpha3<-solve(lambda3*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda4<-.05
alpha4<-solve(lambda4*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})


grd<-seq(0,1,length=100)
K_fit<-outer(grd,t,FUN=myfun)
muhat1<-K_fit%*%alpha1
muhat2<-K_fit%*%alpha2
muhat3<-K_fit%*%alpha3
muhat4<-K_fit%*%alpha4


par(mar=c(2,2,1,1))
plot(c(T),c(Y),xlab="",ylab="")
points(grd,muhat1,type="l",lwd=2,col=1)
points(grd,muhat2,type="l",lwd=2,col=2)
points(grd,muhat3,type="l",lwd=2,col=3)
points(grd,muhat4,type="l",lwd=2,col=4)
leg<-paste("lambda=",c(1,.5,.1,.05))
legend("topright",legend=leg,col=1:4)



y<-c(Y); t<-c(T)
myfun<-function(t,s){
  (1+sqrt(3)*abs(t-s))*exp(-sqrt(3)*abs(t-s))
}
myfun<-Vectorize(myfun)
K<-outer(t,t,FUN=myfun)

  lambda1<-1
alpha1<-solve(lambda1*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda2<-.5
alpha2<-solve(lambda2*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda3<-.1
alpha3<-solve(lambda3*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda4<-.05
alpha4<-solve(lambda4*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})


grd<-seq(0,1,length=100)
K_fit<-outer(grd,t,FUN=myfun)
muhat1<-K_fit%*%alpha1
muhat2<-K_fit%*%alpha2
muhat3<-K_fit%*%alpha3
muhat4<-K_fit%*%alpha4


par(mar=c(2,2,1,1))
plot(c(T),c(Y),xlab="",ylab="")
points(grd,muhat1,type="l",lwd=2,col=1)
points(grd,muhat2,type="l",lwd=2,col=2)
points(grd,muhat3,type="l",lwd=2,col=3)
points(grd,muhat4,type="l",lwd=2,col=4)
leg<-paste("lambda=",c(1,.5,.1,.05))
legend("topright",legend=leg,col=1:4)

y<-c(Y); t<-c(T)
myfun<-function(t,s){
  exp(sin(pi*abs(t-s))^2)
}
myfun<-Vectorize(myfun)
K<-outer(t,t,FUN=myfun)


lambda1<-1
alpha1<-solve(lambda1*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda2<-.5
alpha2<-solve(lambda2*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda3<-.1
alpha3<-solve(lambda3*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})
lambda4<-.05
alpha4<-solve(lambda4*K + t(K)%*%K,t(K)%*%y,tol=10^{-23})


grd<-seq(0,1,length=100)
K_fit<-outer(grd,t,FUN=myfun)
muhat1<-K_fit%*%alpha1
muhat2<-K_fit%*%alpha2
muhat3<-K_fit%*%alpha3
muhat4<-K_fit%*%alpha4


par(mar=c(2,2,1,1))
plot(c(T),c(Y),xlab="",ylab="")
points(grd,muhat1,type="l",lwd=2,col=1)
points(grd,muhat2,type="l",lwd=2,col=2)
points(grd,muhat3,type="l",lwd=2,col=3)
points(grd,muhat4,type="l",lwd=2,col=4)
leg<-paste("lambda=",c(1,.5,.1,.05))
legend("topright",legend=leg,col=1:4)


