##### Lecture 3 Part 1 #####

library(refund); library(ggplot2)
library(dplyr); library(reshape2)
set.seed(2016)

n = 1000;  grid = seq(0, 1, length = 101)
beta1 = sin(grid * 2 * pi)
beta2 = -dnorm(grid, mean=.2, sd=.03) + 
  3*dnorm(grid, mean=.5, sd=.04)+
  dnorm(grid, mean=.75, sd=.05)


par(mfrow=c(1,2),mar=c(2,2,1,1))
plot(grid,beta1,type="l")
plot(grid,beta2,type="l")


X <- matrix(0, nrow=n, ncol=length(grid))
for(i2 in 1:n){
  X[i2,]=X[i2,]+rnorm(length(grid), 0, 1)
  X[i2,]=X[i2,]+runif(1, 0, 5)
  X[i2,]=X[i2,]+rnorm(1, 1, 0.2)*grid
  for(j2 in 1:10){
    e=rnorm(2, 0, 1/j2^(2))
    X[i2,]=X[i2,]+e[1]*sin((2*pi)*grid*j2)
    X[i2,]=X[i2,]+e[2]*cos((2*pi)*grid*j2)
  } }

par(mar=c(2,2,1,1))
matplot(grid,t(X[1:10,]),type="l")

Y1 = X %*% beta1 * .01 + rnorm(n, 0, .4)
Y2 = X %*% beta2 * .01 + rnorm(n, 0, .4)
fit.fpcr1 = pfr(Y1 ~ fpc(X,pve=.85)) 
fit.fpcr2 = pfr(Y2 ~ fpc(X,pve = 0.85))
fit.lin1 = pfr(Y1 ~ lf(X, bs = "ps", k=5, fx=TRUE))
fit.lin2 = pfr(Y2 ~ lf(X, bs = "ps", k=20, fx=TRUE))
# "ps" stands for  "penalized splines", 
# sp = -1 means no penalty is used
fit.pfr1 = pfr(Y1 ~ lf(X, bs = "ps",k=50))
fit.pfr2 = pfr(Y2 ~ lf(X, bs = "ps",k=50))
# if sp is not specified, data 
# driven smoothing is used

par(mar=c(2,2,1,1))
coefs = data.frame(grid = grid,
                   FPCR = coef(fit.fpcr1)$value,
                   Basis = coef(fit.lin1)$value,
                   Penalized = coef(fit.pfr1)$value,
                   Truth = beta1)
# melt stacks the different functions
coefs.m = melt(coefs, id = "grid")
colnames(coefs.m) = c("grid", "Method", "Value")
ggplot(coefs.m, aes(x = grid, y = Value, 
                    color = Method)) + geom_path() + theme_bw()


par(mar=c(2,2,1,1))
coefs = data.frame(grid = grid,
                   FPCR = coef(fit.fpcr1)$value,
                   Basis = coef(fit.lin1)$value,
                   Penalized = coef(fit.pfr1)$value,
                   Truth = beta1)
coefs.m = melt(coefs, id = "grid")
colnames(coefs.m) = c("grid", "Method", "Value")
ggplot(coefs.m, aes(x = grid, y = Value, color = Method)) + geom_path() + theme_bw()

par(mar=c(2,2,1,1))
coefs = data.frame(grid = grid,
                   FPCR = coef(fit.fpcr2)$value,
                   Basis = coef(fit.lin2)$value,
                   Penalized = coef(fit.pfr2)$value,
                   Truth = beta2)
coefs.m = melt(coefs, id = "grid")
colnames(coefs.m) = c("grid", "Method", "Value")
ggplot(coefs.m, aes(x = grid, y = Value, color = Method)) + geom_path() + theme_bw()

###################################
library(fda)
library(refund)

par(mar=c(2,2,1,1))
set.seed(2016)
n = 1000;  grid = seq(0, 1, length = 101)
beta1 = sin(grid * 2 * pi)
beta2 = -dnorm(grid, mean=.2, sd=.03) + 
  3*dnorm(grid, mean=.5, sd=.04)+
  dnorm(grid, mean=.75, sd=.05)
X <- matrix(0, nrow=n, ncol=length(grid))
for(i2 in 1:n){
  X[i2,]=X[i2,]+rnorm(length(grid), 0, 1)
  X[i2,]=X[i2,]+runif(1, 0, 5)
  X[i2,]=X[i2,]+rnorm(1, 1, 0.2)*grid
  for(j2 in 1:10){
    e=rnorm(2, 0, 1/j2^(2))
    X[i2,]=X[i2,]+e[1]*sin((2*pi)*grid*j2)
    X[i2,]=X[i2,]+e[2]*cos((2*pi)*grid*j2)
  } }
Y1 = X %*% beta1 * .01 + rnorm(n, 0, .4)
Y2 = X %*% beta2 * .01 + rnorm(n, 0, .4)
matplot(grid,t(X[1:10,]),type="l")


par(mar=c(2,2,1,1))
mybasis<-create.bspline.basis(c(0,1),nbasis=100)
lambda_all<-10^(seq(0,-4,length=10)); gcv<-numeric(0)
for(lambda in lambda_all){
  myfdpar<-fdPar(mybasis,2,lambda=lambda)
  X.f<-smooth.basis(grid,t(X),myfdpar)
  gcv<-c(gcv,mean(X.f$gcv))}
plot(-log(lambda_all,base=10),gcv)


par(mar=c(2,2,1,1))
lambda=lambda_all[which.min(gcv)]
myfdpar<-fdPar(mybasis,2,lambda=lambda)
X.f<-smooth.basis(grid,t(X),myfdpar)
plot(X.f$fd[2])
points(grid,X[2,])


X.pc<-pca.fd(X.f$fd,3)
cumsum(X.pc$varprop)
lm_fit1<-lm(Y1~X.pc$scores)
lm_fit2<-lm(Y2~X.pc$scores)

print(summary(lm_fit1)$coefficients,digits=2)
anova(lm_fit1)

print(summary(lm_fit2)$coefficients,digits=2)
anova(lm_fit2)

beta_coef1 = coef(X.pc$harmonics)%*%lm_fit1$coef[-1]
beta_hat1<-fd(beta_coef1,mybasis)
beta_coef2 = coef(X.pc$harmonics)%*%lm_fit2$coef[-1]
beta_hat2<-fd(beta_coef2,mybasis)
pfr_fit1<-pfr(Y1~fpc(X,ncomp=3))
pfr_fit2<-pfr(Y2~fpc(X,ncomp=3))


par(mar=c(2,2,1,1))
plot(grid,beta1,type='l')
plot(beta_hat1,add=TRUE,col='red')
points(grid,coef(pfr_fit1)$value,type="l",col="blue")
leg<-c("True","Ours","pfr")
legend("topright",legend=leg,col=c("black","red","blue"),lty=1)


par(mar=c(2,2,1,1))
plot(grid,beta2,type='l')
plot(beta_hat2,add=TRUE,col='red')
points(grid,coef(pfr_fit2)$value,type="l",col="blue")
leg<-c("True","Ours","pfr")
legend("topright",legend=leg,col=c("black","red","blue"),lty=1)

par(mar=c(2,2,1,1))
X.pc<-pca.fd(X.f$fd,6)
lm_fit2<-lm(Y2~X.pc$scores)
beta_coef2 = coef(X.pc$harmonics)%*%lm_fit2$coef[-1]
beta_hat2<-fd(beta_coef2,mybasis)
pfr_fit2<-pfr(Y2~fpc(X,ncomp=6))
plot(grid,beta2,type='l')
plot(beta_hat2,add=TRUE,col='red')
points(grid,coef(pfr_fit2)$value,type="l",col="blue")
leg<-c("True","Ours","pfr")
legend("topright",legend=leg,col=c("black","red","blue"),lty=1)

library(refund)
X = DTI$cca[DTI$case==1,]
Y = DTI$pasat[DTI$case==1]
grid = seq(0,1,length=dim(X)[2])
drop<-unique(which(is.na(X),arr.ind=TRUE)[,1])
X = X[-drop,]; Y = Y[-drop]


k_all = 5:20
aic.fpc = numeric(0)
for(k in k_all){
  pfr.fpc<-pfr(Y~fpc(X,ncomp=k))
  aic.fpc<-c(aic.fpc,pfr.fpc$aic)
}
k = k_all[which.min(aic.fpc)]
pfr.fpc<-pfr(Y~fpc(X,ncomp=k))
pfr.pen<-pfr(Y~lf(X,k=50))
beta_fpc = coef(pfr.fpc)$value
beta_pen = coef(pfr.pen)$value
k


par(mar=c(2,2,1,1))
plot(c(grid),c(beta_fpc),type="l",xlab="",ylab="")
points(grid,beta_pen,type="l",col="blue")
leg = paste(c("FPC, AIC =","Pen, AIC ="),c(round(pfr.fpc$aic),round(pfr.pen$aic)))
legend("bottom",legend=leg,col=c("black","blue"),lty=1)

rbind(summary(pfr.pen)[[22]],summary(pfr.pen)[[24]])

pfr.pen<-pfr(Y~af(X,argvals=grid,k=c(10,10)))
pfr.pen$aic
rbind(summary(pfr.pen)[[22]], summary(pfr.pen)[[24]])

par(mar=c(2,2,1,1))
plot(pfr.pen,scheme=2,xlab="t",ylab="x",main="f(x,t)")


##### Lecture 3 Part 2 #####

library(refund); library(fda)
Y = DTI$cca
X = DTI$case
grid = seq(0,1,length=dim(Y)[2])
drop<-unique(which(is.na(Y),arr.ind=TRUE)[,1])
Y = Y[-drop,]; X = X[-drop]

par(mfrow=c(1,2),mar=c(2,2,1,1))
matplot(grid,t(Y[X==0,]),type="l",col="gray",ylab="",xlab="",ylim=c(.3,.8),main="Control")
points(grid,colMeans(Y[X==0,]),lwd=2,type="l")
matplot(grid,t(Y[X==1,]),type="l",col="gray",ylab="",xlab="",ylim=c(.3,.8),main="Case")
points(grid,colMeans(Y[X==1,]),lwd=2,type="l")

mybasis<-create.bspline.basis(c(0,1),nbasis=20)
Y.f<-Data2fd(grid,t(Y),mybasis)
X.design = cbind(1,X)
fosr.fit<-fosr(fdobj = Y.f, X = X.design, method="OLS")

plot(fosr.fit)


fosr_p <- fosr.perm(fdobj = Y.f, X = X.design, 
                    method="OLS",argvals=grid,nperm=50,prelim=1,plot=FALSE)
fosr_t <- fosr.perm.test(fosr_p)

CI = apply(fosr_t$F.perm,2,quantile,probs=.95)
plot(grid,fosr_t$F,type="l",ylab="",xlab="")
points(grid,CI,type="l",col="red",lty=2)



