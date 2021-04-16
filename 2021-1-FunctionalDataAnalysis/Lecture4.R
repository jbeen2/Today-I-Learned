##### Lecture 4 Part 1 #####
library(fda)
library(refund)
set.seed(2016)

library(RandomFields)

N = 1000
m = 50
pts = seq(0,1,length=m)
RFModel = RMmatern(nu=3/2,var=1,scale=1/4)  # covariance function 
X = RFsimulate(model=RFModel,x=pts,n=N)
X = as.matrix(X)

par(mar=c(3,3,1,1))
matplot(pts,X,type="l")  # X 

RFModel = RMmatern(nu=1/2,var=1,scale=1/4)
Eps = RFsimulate(model=RFModel,x=pts,n=N)
Eps = as.matrix(Eps)

psi_f = function(t,s){5*sin(t*s*pi/2)}  # bivariate 
psi_m = outer(pts,pts,FUN=psi_f)  # bivariate cov function 

library(plot3D)
persp3D(pts,pts,psi_m)  # real beta 

par(mar=c(2,2,1,1))
Y = psi_m%*%X/m + Eps  # (50,50) %*% (50,1000)
matplot(pts,Y,type="l")  # Y 

par(mar=c(2,2,2,1))
Y = t(Y); X = t(X)
pffr_fit = pffr(Y~ff(X,xind=pts),yind=pts)
plot(pffr_fit,select=1)  # pred beta 

par(mar=c(2,2,1,1))
plot(pffr_fit,select=2,pers=TRUE)


# Functional Principal Component Analysis 

# manually...... 
library(fda)
mybasis = create.bspline.basis(c(0,1),nbasis=50)
X.f = Data2fd(pts,t(X),mybasis)
Y.f = Data2fd(pts,t(Y),mybasis)
X.pca = pca.fd(X.f,nharm=3)
Y.pca = pca.fd(Y.f,nharm=5)
psi_hat_m = t(Y.pca$scores)%*%
  (X.pca$scores%*%diag(1/X.pca$values[1:3]))/N   # slide 7p

par(mar=c(2,2,1,1))
psi_hat_coef = coef(Y.pca$harmonics)%*%
  psi_hat_m%*%t(coef(X.pca$harmonics))
psi_hat_f = bifd(psi_hat_coef,mybasis,mybasis)
psi_hat_ev = eval.bifd(pts,pts,psi_hat_f)
persp3D(pts,pts,psi_hat_ev)

Delta = t(Y.pca$scores)%*%(X.pca$scores)/N
norm = Y.pca$values[1:5]%*%t(X.pca$values[1:3])
T_stat = N*sum(Delta^2/norm)
T_stat
pchisq(T_stat,df = 5*3,lower.tail=FALSE)

###################################
# Functional GLM
###################################

##### Lecture 4 Part 2 #####

library(fda)
library(refund)
set.seed(2016)

library(fda)
N<-1000; M<-50
pts = seq(0,1,length=M)
alpha = 1
c = 1
beta = -c*sin(2*pi*pts) 

library(RandomFields)
RFModel = RMmatern(nu=5/2,var=1,scale=1/4)
X = RFsimulate(model=RFModel,x=pts,n=N)
X = as.matrix(X)
X = X + pts*3
X = t(X)

par(mar=c(2,2,0,0))
matplot(pts,t(X),type="l")

Z = (alpha + X%*%beta/M + rnorm(N,0,1))
Y = Z>0

k_all = 4:10
aic_all = numeric(0)
for(k in k_all){
  pfr.fit<-pfr(Y~lf(X,k=k,fx=TRUE),
               family=binomial(link="probit"))
  tmp = extractAIC(pfr.fit)[2]
  aic_all = c(aic_all,tmp)
}


par(mar=c(2,2,0,0))
k = k_all[which.min(aic_all)]
pfr.fit<-pfr(Y~lf(X,k=k,fx=TRUE),
             family=binomial(link="probit"))
plot(pfr.fit)
points(pts,beta,type="l",col='red')

par(mar=c(2,2,0,0))
library(refund)
pfr.fit<-pfr(Y~lf(X,k=50),
             family=binomial(link="probit"))
plot(pfr.fit)
points(pts,beta,type="l",col='red')




library(fda); library(RandomFields)
N<-1000; M<-50
pts = seq(0,1,length=M)
pts_grid = outer(pts,pts,FUN="*")
alpha = 0; c = 1
beta = c*sin(pi*pts_grid) 
RFModel = RMmatern(nu=5/2,var=1,scale=1/4)
X = RFsimulate(model=RFModel,x=pts,n=N)
Eps = RFsimulate(model=RFModel,x=pts,n=N)
X = as.matrix(X)
Eps = .1* as.matrix(Eps); Eps = t(Eps)
X = X + pts*3 - 3/2; X = t(X)
Z = alpha + X%*%beta/M + Eps
Y = Z>0


par(mfrow=c(2,2),mar=c(2,2,1,1))
for(i in 1:4){
  plot(pts,Y[i,],ylim=c(0,1))
}

pffr.fit<-pffr(Y~ff(X,xind=pts),
               family=binomial(link="probit"),yind=pts)

par(mar=c(2,2,1,1),mfrow=c(1,2))
plot(pffr.fit,select=1)
plot(pffr.fit,select=2,pers=TRUE)


bpars = list(bs="ps",k= 4,fx=TRUE)
intpar = list(bs = "ps", k = 5, fx=TRUE)
pffr.fit2<-pffr(Y~ff(X,xind=pts,splinepars=bpars),
                bs.int=intpar, family=binomial(link="probit"),yind=pts)


par(mar=c(2,2,1,1))
plot(pffr.fit2,select=2,pers=TRUE)







