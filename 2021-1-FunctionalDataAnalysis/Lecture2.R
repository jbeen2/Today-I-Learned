##### Lecture 2 Part 1 #####

library(fda)
library(ggplot2)

set.seed(201609); library(fields); library(expm)
m<-100; times<-seq(0,1,length=m)
range<-1; nu1=1/2; sig2<-1; nu2=3/2
Matern(.5,range=range,nu=nu1)

# Gaussian Process Matern Function 
d_mat<-abs(outer(times,times,"-"))  # Gaussian Process, matrix 만들기 
C_1<-apply(d_mat,c(1,2),FUN=Matern,range=range,nu=nu1)
C_1<-C_1*sig2
C_1_sq<-sqrtm(C_1)
C_2<-apply(d_mat,c(1,2),FUN=Matern,range=range,nu=nu2)
C_2<-C_2*sig2
C_2_sq<-sqrtm(C_2)
Z<-rnorm(m)
X1<-C_1_sq%*%Z; X2<-C_2_sq%*%Z

par(mar=c(3,3,1,1),mfrow=c(1,2))
plot(times,X1,type="l")
plot(times,X2,type="l")

# manually ver... 
par(mar=c(3,3,1,1),mfrow=c(1,2))
Xd1<-(tail(X1,-1) - head(X1,-1))*m  # 기울기, 미분값 만들기 
Xd2<-(tail(X2,-1) - head(X2,-1))*m
plot(times[-1],Xd1,type="l")
plot(times[-1],Xd2,type="l")


mybasis<-create.bspline.basis(c(0,1),nbasis=50)
X2.f<-Data2fd(times,X2,mybasis)
X2d.f<-deriv.fd(X2.f)  # derivative 

par(mar=c(3,3,1,1),mfrow=c(2,2))
plot(times,X2,type="l",main="Numeric");
plot(X2.f,main="Basis");
plot(times[-1],Xd2,type="l");
plot(X2d.f)


times = growth$age; GHeight = growth$hgtf
my_basis<-create.bspline.basis(c(1,18),
                               nbasis=10,norder=5)
my_par<-fdPar(my_basis,Lfdobj=2,lambda=1)
GHeight.S<-smooth.basis(times,GHeight,my_par)
names(GHeight.S)[1:5]; names(GHeight.S)[6:9]



par(mar=c(3,3,1,1))
mean(GHeight.S$gcv)
hist(GHeight.S$gcv)


par(mar=c(3,3,1,1),mfrow=c(2,2))

main1=paste("lambda=",1," GCV=",round(mean(GHeight.S$gcv),3))
plot(GHeight.S$fd,main=main1)

lambda=0.05
my_par<-fdPar(my_basis,Lfdobj=2,lambda=lambda)  # penalized 된 것  
GHeight.S2<-smooth.basis(times,GHeight,my_par)
main2=paste("lambda=",lambda," GCV=",round(mean(GHeight.S2$gcv),3))
plot(GHeight.S2$fd,main=main2)

lambda=0.01
my_par<-fdPar(my_basis,Lfdobj=2,lambda=lambda)
GHeight.S3<-smooth.basis(times,GHeight,my_par)
main3=paste("lambda=",lambda," GCV=",round(mean(GHeight.S3$gcv),3))
plot(GHeight.S3$fd,main=main3)

lambda=20000
my_par<-fdPar(my_basis,Lfdobj=2,lambda=lambda)
GHeight.S4<-smooth.basis(times,GHeight,my_par)
main4=paste("lambda=",lambda," GCV=",round(mean(GHeight.S4$gcv),3))
plot(GHeight.S4$fd,main=main4)


par(mar=c(3,3,1,1),mfrow=c(2,2))
GHeight.F<-GHeight.S$fd
GHeightD.F<-deriv.fd(GHeight.F,Lfdobj=2)
muD.F<-mean.fd(GHeightD.F)
plot(GHeightD.F,col='grey',main=main1);plot(muD.F,lwd=2,add=TRUE)

GHeight.F<-GHeight.S2$fd
GHeightD.F<-deriv.fd(GHeight.F,Lfdobj=2)
muD.F<-mean.fd(GHeightD.F)
plot(GHeightD.F,col='grey',main=main2);plot(muD.F,lwd=2,add=TRUE)

GHeight.F<-GHeight.S3$fd
GHeightD.F<-deriv.fd(GHeight.F,Lfdobj=2)
muD.F<-mean.fd(GHeightD.F)
plot(GHeightD.F,col='grey',main=main3);plot(muD.F,lwd=2,add=TRUE)

GHeight.F<-GHeight.S4$fd
GHeightD.F<-deriv.fd(GHeight.F,Lfdobj=2)
muD.F<-mean.fd(GHeightD.F)
plot(GHeightD.F,col='grey',main=main4);plot(muD.F,lwd=2,add=TRUE)


nbasis = 365; yearRng = c(0,365)
daybasis = create.fourier.basis(yearRng, nbasis)
logprecav = CanadianWeather$dailyAv[,,'log10precip']
dayprecfd <- with(CanadianWeather, smooth.basis(day.5,
                                                logprecav, daybasis,
                                                fdnames=list("Day","Station","log10(mm)"))$fd)
names(CanadianWeather)[1:3];names(CanadianWeather)[4:6];names(CanadianWeather)[7:8]


par(mar=c(3,3,1,1))
plot(dayprecfd[1]);points(day.5,logprecav[,1])

par(mar=c(3,3,1,1))
daybasis = create.fourier.basis(yearRng, nbasis)
mypar = fdPar(daybasis,Lfdobj=2,lambda=10000)
Y.f<-smooth.basis(day.5,logprecav,mypar)
plot(dayprecfd[1],col='grey');points(day.5,logprecav[,1])
plot(Y.f$fd[1],add=TRUE,lwd=2)

mypar = fdPar(daybasis,Lfdobj=2,lambda=.001)
Y.f2<-smooth.basis(day.5,logprecav,mypar)
plot(Y.f2$fd[1],add=TRUE,lwd=2)
mean(Y.f2$gcv)
mean(Y.f$gcv)

##### Lecture 2 Part 2 #####

# sin 함수 : 갖고 싶은 올바른 형태 갖지 못함 
par(mar=c(3,3,1,1))
myfun<-function(a,x){
  if(2*pi*x - a < 0){return(0)
  }else if(2*pi*x - a > 2*pi){return(0)
  }else{return(sin(2*pi*x-a))}
}
times<-seq(0,1.6,length=100)
x1<-sapply(times,FUN=myfun,a=3/4)
x2<-sapply(times,FUN=myfun,a=6/4)
x3<-sapply(times,FUN=myfun,a=9/4)
x4<-sapply(times,FUN=myfun,a=12/4)
mu<-(x1+x2+x3+x4)/4
plot(times,x1,type="l",ylab="",xlab="")
points(times,x2,type="l")
points(times,x3,type="l")
points(times,x4,type="l")
points(times,mu,type="l",lty=2)


par(mar=c(3,3,1,1))
plot(times,x1,type="l",ylab="",xlab="")
num_loc<-2  # 내가 원하는 landmark click ! 
t_loc<-locator(num_loc)  # 내가 찍은 landmark 의 좌표가 나옴 


# curve 를 돌면서 align  
hgtbasis <- create.bspline.basis(c(1,18), norder = 6, breaks=growth$age)
growfdPar <- fdPar(hgtbasis, 4, 10^(-.5))
X.f<-smooth.basis(growth$age,growth$hgtm,growfdPar)$fd
X.reg<-register.fd(X.f)  # register.fd 함부로 하면 위험한 이유...?? 

names(X.reg)


par(mar=c(3,3,1,1));par(mfrow=c(1,2))
plot(X.f,col='gray'); plot(mean(X.f),lwd=2,add=TRUE)
plot(X.reg$regfd,col='gray'); plot(mean(X.reg$regfd),lwd=2,add=TRUE)


X.fd<-deriv(X.f,2)
X.fd2<-register.fd(X.fd)$regfd


# alignment 전 vs 후 
par(mar=c(3,3,1,1));par(mfrow=c(1,2))
plot(X.fd,col='gray',ylim=c(-4,2)); plot(mean(X.fd),lwd=2,add=TRUE)
plot(X.fd2,col='gray',ylim=c(-4,2)); plot(mean(X.fd2),lwd=2,add=TRUE)


library(refund)
Corp<-DTI$cca
drop<-unique(which(is.na(Corp),arr.ind=TRUE)[,1])
Corp<-Corp[-drop,] # Missing value
pts<-seq(0,1,length=93)
my_basis<-create.bspline.basis(c(0,1),
                               nbasis=100,norder=6)

lambda_all<-10^(-(10:20)/2)  # lambda 의 candidate 를 만들어 줌 
gcv_all<-numeric(0)
for(lambda in lambda_all){
  myPar<-fdPar(my_basis,2,lambda)  # i 번째 candidate 로 fdPar 
  Corp.F<-smooth.basis(pts,t(Corp),myPar)
  gcv_all<-c(gcv_all,mean(Corp.F$gcv))}  # mean 을 가장 작게 만들어 주는 lambda 값 찾기 

# 최적의 lambda 찾기 
par(mar=c(5,3,1,1))
lambda_all[which.min(gcv_all)]  
plot(-log10(lambda_all),gcv_all)

best_lam<-lambda_all[which.min(gcv_all)]  # 제일 좋은 lambda 로 
myPar<-fdPar(my_basis,2,best_lam)
Corp.F<-smooth.basis(pts,t(Corp),myPar)  # smoothing 
plot(Corp.F,main=paste("lambda=",best_lam))


par(mar=c(3,3,1,1))
best_lam<-.00001  # lambda=1000000 : 직선의 형태 
myPar<-fdPar(my_basis,2,best_lam)
Corp.F<-smooth.basis(pts,t(Corp),myPar)$fd
plot(Corp.F,main=paste("lambda=",best_lam))

Corp.FR<-register.fd(Corp.F)$regfd  # Functional data 에 대해 registeration  
Corp.FD<-deriv(Corp.F)
Corp.FDR<-register.fd(Corp.FD)$regfd  # 미분값에 대해 registration 

# Unaligned vs Aligned 
par(mar=c(3,3,1,1),mfrow=c(1,2))
plot(Corp.F,main="Unaligned",col='gray');plot(mean.fd(Corp.F),add=TRUE,lwd=3)
plot(Corp.FR,main="Aligned",col='gray');plot(mean.fd(Corp.FR),add=TRUE,lwd=3)

par(mar=c(3,3,1,1),mfrow=c(1,2))
plot(Corp.FD,main="Unaligned",col='gray');plot(mean.fd(Corp.FD),add=TRUE,lwd=3)
plot(Corp.FDR,main="Aligned",col='gray');plot(mean.fd(Corp.FDR),add=TRUE,lwd=3)




