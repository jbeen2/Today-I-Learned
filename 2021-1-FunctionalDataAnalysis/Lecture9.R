##### Lecture 9 Part 1 #####
library(fda)
library(refund)
set.seed(2016)

# VAS Longitudinal Measurements
Data<-read.csv("va.csv")
# VAS Baseline Measurements
Base<-read.csv("bv.csv")
# Other Subject Information
general_info<-read.csv("gi.csv")

# Each row corresponds to a particular subject,
# so we start by pulling off the subject 'ids' for each dataset.
IdData<-Data[,1]
IdBase<-Base[,1]
IdGeneral<-general_info[,1]
# Unique ids
Code_levels<-levels(factor(IdBase))
# Number of subects : 1185 -> 1185^(1/4)=5.8, 어떤 데이터는 sparse (6개), 어떤 데이터는 non sparse 
n<-length(Code_levels)   # -> Sparse Functional Data 로 보겠다 


# The "time" variable which is given by the week of the study.
Week<-Data[,"week"]
# Unique weeks (=unique week, 4 ~ 104)
Week_levels<-levels(factor(Week))
# Select Visual Acuity values
VA<-Data[,"studyeye_va"]
VAB<-Base[,"studyeye_va"]  # 기준점 score 

Age<-general_info[,"age"]
Y_l<-by(c(VAB,VA),c(IdBase,IdData),c) # 각 id 에 대해 관측된 점수들 
T_l<-by(c(rep(0,times=n),Week),c(IdBase,IdData),c)
counts<-sapply(Y_l,length,simplify=TRUE)
# Drop those with only a baseline measurement
Y_l<-Y_l[counts!=1]; T_l<-T_l[counts!=1]   # baseline 만 있는 경우 제외  
counts<-counts[counts!=1]
n<-length(counts)  # 1177 명의 실험 참가자! 

par(mar=c(2,2,1,1))
hist(counts)

par(mar=c(2,2,1,1))
Y<-unlist(Y_l); T<-unlist(T_l)
plot(T,Y)  # Y : 점수, T : 주 

par(mar=c(2,2,1,1))
loess_fit<-loess(Y~T,span = 1, degree=2)  # loess : local polynomial regression (Lecture8.30p)
grd<-seq(0,104,length=50)
yhat<-predict(loess_fit,newdata= grd)
plot(T,Y,xlab="Week",ylab="VAS")
points(grd,yhat,type="l",col="red",lwd=3)  # 점수 우상향 

par(mar=c(2,2,1,1))
lmfit<-lm(Y~as.factor(T)-1)
yhat0<-coef(lmfit)
plot(c(0,Week_levels),yhat0)

yhat1<-yhat
loess_fit<-loess(Y~T,span = .75, degree=2)  # span = smoothness 
yhat2<-predict(loess_fit,newdata= grd)
loess_fit<-loess(Y~T,span = .5, degree=2)
yhat3<-predict(loess_fit,newdata= grd)
loess_fit<-loess(Y~T,span = .25, degree=2)  # span parameter 작으면 wiggly (Overfitting)
yhat4<-predict(loess_fit,newdata= grd) 


points(grd,yhat1,type="l",col=2,lwd=3)
points(grd,yhat2,type="l",col=3,lwd=3)
points(grd,yhat3,type="l",col=4,lwd=3)
points(grd,yhat4,type="l",col=5,lwd=3)

leg<-paste("span=",c(1,.75,.5,.25))
legend("bottomright",legend=leg,col=2:5,lty=1)


# GAM 
par(mar=c(2,2,1,1))
library(mgcv)
gam_fit<-gam(Y~s(T,sp=.1))  # gam : generalize additive model (Y = f1(X1)+f2(X2)+f3(X3))
pred_gam<-predict(gam_fit,newdata=list(T=grd))
plot(c(0,Week_levels),yhat0)
points(grd,pred_gam,col=2,lwd=3,type="l")

par(mar=c(2,2,1,1))
gam_fit<-gam(Y~s(T,sp=1))
pred_gam2<-predict(gam_fit,newdata=list(T=grd))

gam_fit<-gam(Y~s(T,sp=10))
pred_gam3<-predict(gam_fit,newdata=list(T=grd))

gam_fit<-gam(Y~s(T,sp=100))
pred_gam4<-predict(gam_fit,newdata=list(T=grd))

plot(c(0,Week_levels),yhat0)
points(grd,pred_gam,col=2,lwd=3,type="l")
points(grd,pred_gam2,col=3,lwd=3,type="l")
points(grd,pred_gam3,col=4,lwd=3,type="l")
points(grd,pred_gam4,col=5,lwd=3,type="l")

leg<-paste("sp=",c(.1,1,10,100))
legend("bottomright",legend=leg,col=2:5,lty=1)


# Covariance Estimation 
muhat<-gam(Y~s(T,sp=1))
Z_l<-vector("list",n)
TT_l<-vector("list",n)
for(i in 1:n){
  tmp<-Y_l[[i]]   # 첫번째 사람의 점수 
  tmpT<-T_l[[i]]  # 첫번째 사람의 timepoint 
  tmp<-tmp - predict(muhat, newdata=
                       data.frame(T=tmpT))  # y - muhat 
  tmp2<-outer(tmp,tmp,"*")
  Z_l[[i]]<-tmp2[lower.tri(tmp2,diag=FALSE)]
  tmpT2<-expand.grid(tmpT,tmpT)
  tmpT3<-tmpT2[c(lower.tri(tmp2,diag=FALSE)),]
  TT_l[[i]]<-tmpT3
}

mydata<-cbind(do.call(rbind,TT_l),unlist(Z_l)) # timepoint 이용해서 smoothing 
z<-mydata[,3] ; x1<-mydata[,1] ; x2<-mydata[,2]
GAMF_C<-gam(z~s(x1,x2,sp = 1000))
m<-50
pts<-seq(0,104,length=m)
Sigma<-matrix(nrow=m,ncol=m)
for(i in 1:m){
  for(j in 1:i){
    Sigma[i,j]<-predict(GAMF_C,newdata =
                          data.frame(x1=pts[i],x2=pts[j]))
    Sigma[j,i]<-Sigma[i,j]
  }}

par(mar=c(2,2,1,1))
library(plot3D)
persp3D(pts,pts,Sigma)

# (BLUP : Best Linear Unbiased Prediction)


##### Lecture 9 Part 2 #####
# PACE : Principal Analysis by Conditional Expectation 
library(fda)
library(refund)
set.seed(2016)

set.seed(2016); library(MASS); library(fdapace)
M_all<-100; M<-5; N<-100; sig2<-.5; rho<-.1;
grid_all<-seq(0,1,length= M_all)

mu <-function(t){
  .5*dnorm(t,mean = 1/3,sd=1/18)+
    .1*dnorm(t,mean = 2/3,sd=1/18)}
C<-function(d){exp(-abs(d)^2/rho)}

# data generate 
mu_all<-mu(grid_all)
Sig_all<-C(outer(grid_all,grid_all,FUN="-")) # outer product 
Y_all<-mvrnorm(N,mu_all,Sig_all) # Xt + subject error 
Y_l<-vector("list",length=N)
T_l<-vector("list",length=N)
for(i in 1:N){
  spts<-sort(sample(1:M_all,size=M,replace=FALSE))
  Y_l[[i]] <-Y_all[i,spts] + rnorm(M,0,sqrt(sig2))  # + measurement error 
  T_l[[i]] <-grid_all[spts]}  # 100^(1/4) = 3 -> Sparse Dataset 


# FPCA 
myopts<-list(dataType="Sparse",nRegGrid=M_all)  # Sparse DataType, Grid = #100 
FPCA.fit<-FPCA(Y_l,T_l,optns=myopts)
names(FPCA.fit)

par(mar=c(2,2,1,1))
muhat<-FPCA.fit$mu; pts<-FPCA.fit$workGrid
plot(unlist(T_l),unlist(Y_l))
points(grid_all,mu_all,type="l",lwd=3,col=2)
points(pts,muhat,type="l",col=5,lwd=3)

library(plot3D)
par(mfrow=c(1,2))
par(mar=c(2,2,1,1))
Sigmahat<-FPCA.fit$smoothedCov
persp3D(pts,pts,Sigmahat,colkey = FALSE)  # Estimate 
persp3D(grid_all,grid_all,Sig_all,colkey = FALSE)  # True 

E<-eigen(Sig_all)
Ehat<-eigen(Sigmahat)  # eigen decomposition 
v1<-E$vectors[,1]*sqrt(length(grid_all))
v1hat<-Ehat$vectors[,1]*sqrt(length(pts))
sgn<-sign(t(v1)%*%v1hat)
v1hat<-sgn*v1hat

# 정확도 비교 : First Principal Component 
par(mar=c(2,2,1,1))
plot(grid_all,v1,type="l",ylim=c(-.5,1.5),main="First E. Function")
points(pts,v1hat,type="l",col=2)

# Second Principal Component 
par(mar=c(2,2,1,1))
v2<-E$vectors[,2]*sqrt(length(grid_all))
v2hat<-Ehat$vectors[,2]*sqrt(length(pts))
sgn<-sign(t(v2)%*%v2hat)
v2hat<-sgn*v2hat
par(mar=c(2,2,1,1))
plot(grid_all,v2,type="l",ylim=c(-1.5,2.5),main="Second E. Function")
points(pts,v2hat,type="l",col=2)


# Sample size (N), timepoints(M) 늘어날수록 더 정확해진다 

# 1. 
# data generate 
M_all<-100; M<-5; N<-200; sig2<-.5; rho<-.1;  # 200^(1/4) = 3.7 
grid_all<-seq(0,1,length= M_all)
mu <-function(t){
  .5*dnorm(t,mean = 1/3,sd=1/18)+
    .1*dnorm(t,mean = 2/3,sd=1/18)}
C<-function(d){exp(-abs(d)^2/rho)}
mu_all<-mu(grid_all)  # True mu 
Sig_all<-C(outer(grid_all,grid_all,FUN="-"))
Y_all<-mvrnorm(N,mu_all,Sig_all)
Y_l<-vector("list",length=N)
T_l<-vector("list",length=N)
for(i in 1:N){
  spts<-sort(sample(1:M_all,size=M,replace=FALSE))
  Y_l[[i]] <-Y_all[i,spts] + rnorm(M,0,sqrt(sig2))
  T_l[[i]] <-grid_all[spts]}

# FPCA 
FPCA.fit<-FPCA(Y_l,T_l,optns=myopts)
par(mar=c(2,2,1,1))
muhat<-FPCA.fit$mu; pts<-FPCA.fit$workGrid  # Estimate mu 

err<-mean((muhat-mu_all)^2) # 0.1 
plot(unlist(T_l),unlist(Y_l),ylim=c(-2,4),main=paste("L2 Error=",round(err,digits=3)))
points(grid_all,mu_all,type="l",lwd=3,col=2)
points(pts,muhat,type="l",col=5,lwd=3)

# Covariance 
par(mfrow=c(1,2))
par(mar=c(2,2,1,1))
Sigmahat<-FPCA.fit$smoothedCov
err<-mean((Sigmahat-Sig_all)^2)  # 0.033648 
persp3D(pts,pts,Sigmahat,colkey = FALSE,zlim=c(0,1.1))
persp3D(grid_all,grid_all,Sig_all,colkey = FALSE,zlim=c(0,1.1))
mtext(paste("L2 Error=",round(err,digits=3)),outer=TRUE,line=-2,cex=1.75)


# 2. 
M_all<-100; M<-5; N<-1000; sig2<-.5; rho<-.1; # obs 올림 
grid_all<-seq(0,1,length= M_all)
mu <-function(t){
  .5*dnorm(t,mean = 1/3,sd=1/18)+
    .1*dnorm(t,mean = 2/3,sd=1/18)}
C<-function(d){exp(-abs(d)^2/rho)}
mu_all<-mu(grid_all)
Sig_all<-C(outer(grid_all,grid_all,FUN="-"))
Y_all<-mvrnorm(N,mu_all,Sig_all)
Y_l<-vector("list",length=N)
T_l<-vector("list",length=N)
for(i in 1:N){
  spts<-sort(sample(1:M_all,size=M,replace=FALSE))
  Y_l[[i]] <-Y_all[i,spts] + rnorm(M,0,sqrt(sig2))
  T_l[[i]] <-grid_all[spts]}

FPCA.fit<-FPCA(Y_l,T_l,optns=myopts)
par(mar=c(2,2,1,1))
muhat<-FPCA.fit$mu; pts<-FPCA.fit$workGrid
err<-mean((muhat-mu_all)^2)  
plot(unlist(T_l),unlist(Y_l),ylim=c(-2,4),main=paste("L2 Error=",round(err,digits=3)))
points(grid_all,mu_all,type="l",lwd=3,col=2)
points(pts,muhat,type="l",col=5,lwd=3)  # L2error = 0.08 : sample size 늘려서 error 작아짐 

par(mfrow=c(1,2))
par(mar=c(2,2,1,1))
Sigmahat<-FPCA.fit$smoothedCov
err<-mean((Sigmahat-Sig_all)^2)
persp3D(pts,pts,Sigmahat,colkey = FALSE,zlim=c(0,1.1))
persp3D(grid_all,grid_all,Sig_all,colkey = FALSE,zlim=c(0,1.1))
mtext(paste("L2 Error=",round(err,digits=3)),outer=TRUE,line=-2,cex=1.75)  # 0.009 


# 3. 
M_all<-100; M<-10; N<-100; sig2<-.5; rho<-.1;  # 100 timepoint, 100 obs -> 더 정확하게 추정 
grid_all<-seq(0,1,length= M_all)
mu <-function(t){
  .5*dnorm(t,mean = 1/3,sd=1/18)+
    .1*dnorm(t,mean = 2/3,sd=1/18)}
C<-function(d){exp(-abs(d)^2/rho)}
mu_all<-mu(grid_all)
Sig_all<-C(outer(grid_all,grid_all,FUN="-"))
Y_all<-mvrnorm(N,mu_all,Sig_all)
Y_l<-vector("list",length=N)
T_l<-vector("list",length=N)
for(i in 1:N){
  spts<-sort(sample(1:M_all,size=M,replace=FALSE))
  Y_l[[i]] <-Y_all[i,spts] + rnorm(M,0,sqrt(sig2))
  T_l[[i]] <-grid_all[spts]}
FPCA.fit<-FPCA(Y_l,T_l,optns=myopts)
par(mar=c(2,2,1,1))
muhat<-FPCA.fit$mu; pts<-FPCA.fit$workGrid
err<-mean((muhat-mu_all)^2)
plot(unlist(T_l),unlist(Y_l),main=paste("L2 Error=",round(err,digits=3)))
points(grid_all,mu_all,type="l",lwd=3,col=2)
points(pts,muhat,type="l",col=5,lwd=3)  # 0.016 

par(mfrow=c(1,2))
par(mar=c(2,2,1,1))
Sigmahat<-FPCA.fit$smoothedCov
persp3D(pts,pts,Sigmahat,colkey = FALSE)
persp3D(grid_all,grid_all,Sig_all,colkey = FALSE)
mtext(paste("L2 Error=",round(err,digits=3)),outer=TRUE,line=-2,cex=1.75)  # 0.106 


# 4. 
M_all<-100; M<-100; N<-100; sig2<-.5; rho<-.1;  # timepoint 올림! = 더 이상 sparse 하지 않음 
grid_all<-seq(0,1,length= M_all)
mu <-function(t){
  .5*dnorm(t,mean = 1/3,sd=1/18)+
    .1*dnorm(t,mean = 2/3,sd=1/18)}
C<-function(d){exp(-abs(d)^2/rho)}
mu_all<-mu(grid_all)
Sig_all<-C(outer(grid_all,grid_all,FUN="-"))
Y_all<-mvrnorm(N,mu_all,Sig_all)
Y_l<-vector("list",length=N)
T_l<-vector("list",length=N)
for(i in 1:N){
  spts<-sort(sample(1:M_all,size=M,replace=FALSE))
  Y_l[[i]] <-Y_all[i,spts] + rnorm(M,0,sqrt(sig2))
  T_l[[i]] <-grid_all[spts]}
FPCA.fit<-FPCA(Y_l,T_l,optns=myopts)
par(mar=c(2,2,1,1))
muhat<-FPCA.fit$mu; pts<-FPCA.fit$workGrid
err<-mean((muhat-mu_all)^2)
plot(unlist(T_l),unlist(Y_l),main=paste("L2 Error=",round(err,digits=3)))
points(grid_all,mu_all,type="l",lwd=3,col=2)
points(pts,muhat,type="l",col=5,lwd=3)

par(mfrow=c(1,2))
par(mar=c(2,2,1,1))
Sigmahat<-FPCA.fit$smoothedCov
persp3D(pts,pts,Sigmahat,colkey = FALSE)
persp3D(grid_all,grid_all,Sig_all,colkey = FALSE)
mtext(paste("L2 Error=",round(err,digits=3)),outer=TRUE,line=-2,cex=1.75)



# VAS Longitudinal Measurements
Data<-read.csv("data tables/va.csv")
# VAS Baseline Measurements
Base<-read.csv("data tables/bv.csv")
# Other Subject Information
general_info<-read.csv("data tables/gi.csv")
# Each row corresponds to a particular subject,
# so we start by pulling off the subject ids
# for each dataset.
IdData<-Data[,1]
IdBase<-Base[,1]
IdGeneral<-general_info[,1]
# Unique ids
Code_levels<-levels(factor(IdBase))
# Number of subects
n<-length(Code_levels)

# The "time" variable which is given
# by the week of the study.
Week<-Data[,"week"]
# Unique weeks
Week_levels<-levels(factor(Week))
# Select Visual Acuity values
VA<-Data[,"studyeye_va"]
VAB<-Base[,"studyeye_va"]

Age<-general_info[,"age"]
Y_l<-by(c(VAB,VA),c(IdBase,IdData),c)
T_l<-by(c(rep(0,times=n),Week),c(IdBase,IdData),c)
counts<-sapply(Y_l,length,simplify=TRUE)
# Drop those with only a baseline measurement
Y_l<-Y_l[counts!=1]; T_l<-T_l[counts!=1]
Age<-Age[counts!=1]
counts<-counts[counts!=1]
n<-length(counts)


# 감으로 잡지 않고 바로 사용 가능! 
library(fdapace) 
FPCA_fit<-FPCA(Y_l,T_l,optns=list(dataType="Sparse")) # but 겁나 오래걸림.. 
muhat<-FPCA_fit$mu; pts<-FPCA_fit$workGrid

par(mar=c(2,2,1,1))
plot(unlist(T_l),unlist(Y_l),ylim=c(40,80))
points(pts,muhat,type="l",col="red",lwd=3)


v_hat<-FPCA_fit$phi
lambda_hat<-FPCA_fit$lambda
scores<-FPCA_fit$xiEst
pts<-FPCA_fit$workGrid
lm_fit<-lm(scores~Age)
alpha_hat<-lm_fit$coefficients[1,]
beta_hat<-lm_fit$coefficients[2,]
beta_hat_fun<-v_hat%*%beta_hat

par(mar=c(2,2,1,1))
plot(pts,beta_hat_fun,type="l")


Y<-unlist(Y_l)
T<-unlist(T_l)
subj<-rep(1:n,times=counts)
Age_tmp = rep(Age,times=counts)
Age_tmp<-Age_tmp[!is.na(Y)]
T<-T[!is.na(Y)]; subj<-subj[!is.na(Y)]; Y<-Y[!is.na(Y)]
Y_all<-data.frame(.obs=subj,.index=T,.value=Y)
Xdata=data.frame(X0=rep(1,times=n),X=Age)
pffr_fit<-pffr(Y~X0+X, data=Xdata, ydata=Y_all,
               bs.yindex = list(bs = "ps",k = 20, m = c(2, 1)))

par(mar=c(2,2,1,1))
plot(pffr_fit,select=2,ylim=c(-2,2))
