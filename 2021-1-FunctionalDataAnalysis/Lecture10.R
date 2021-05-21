library(RCurl);
library(demography);
library(MortalitySmooth);
library(ftsa);
library(vars)


usa <- hmd.mx("USA", "kim.kyongwon@yahoo.com", "1234", "USA")
usa1950 <- extract.years(usa, years=1950:2010)

smus <- smooth.demogdata(usa1950)
plot(usa1950, years=2003, series="male", type="p", pch=1)
lines(smus, years=2003, series="male")


fdm.male <- fdm(smus, series="male", order=3)
forecast.fdm.male <- forecast(fdm.male, h=30)
plot(forecast.fdm.male, plot.type="component")
plot(forecast.fdm.male)


############################################################################
# Another example
############################################################################

library("ftsa"); par(mfrow=c(2,1))
# Plot the historical data in gray
plot(Australiasmoothfertility, col = gray(0.8), xlab = "Age", ylab = "Births per 1,000 females",main = "Forecasted fertility rates (2007-2026)")
# Plot the forecasts in rainbow color
plot(forecast(ftsm(Australiasmoothfertility, order = 2), h =20), add = TRUE)
legend("topright", c("2007", "2026"), col = c("red", "blue"), lty = 1)
# Repeat the above steps but use recursive forecasts
plot(Australiasmoothfertility,col=gray(0.8),xlab="Age",ylab = "Births per 1,000 females",main = "Recursive Forecasts")
plot(ftsmiterativeforecasts(Australiasmoothfertility, components = 2, iteration = 20), add = TRUE)
legend("topright", c("2007", "2026"), col = c("red", "blue"), lty = 1)





###########################################################
# Generation and estimation of the FAR(1) model
###########################################################


library(fda)
m=100 # each function is observed at m+1 points, including 0 and 1
burnin=50 # the first 50 functions are a burn in period
N=200 # number of functions to simulate
N1=N+burnin
alpha=9/4
# Create 2 matrices, whose entries are 0s.
# Each column represents one function.
X<- matrix(rep(0, (m+1)*N1),m+1,N1)

epsilon<- matrix(rep(0, (m+1)*N1),m+1,N1)
epsilon[,1]<-rnorm(1)*sin(pi*(0:m/m))+0.5*rnorm(1)*cos(2*pi*(0:m/m))
# the following loop simulates FAR(1).
for(i in 2:N1){
  epsilon[,i]<-rnorm(1)*sin(pi*(0:m/m))+0.5*rnorm(1)*cos(2*pi*(0:m/m))
  X[,i]<-alpha*(1/m)^2*sum((1:m)*X[2:(m+1),i-1])*(0:m/m)+epsilon[,i]
}
X=X[,-(1:burnin)] # Remove the burn in period functions


last=100
plot.ts(c(X[,(N-last+1):N]),ylim=c(min(X[,(N-last):N])-0.5,0.5 +max(X[,(N-last):N])),axes=F,xlab="",ylab="",lwd=2)
axis(2); axis(1,tick=F,labels=F); abline(h=0)
abline(v=seq(0,last*(m+1),by=m+1), lty=2); box()




basisfd=10 # number of basis functions to represent each functional observation
basis=create.bspline.basis(c(0,1),nbasis=basisfd,norder=4)
fdX=Data2fd(argvals=0:m/m,X, basis)


p=4 # number of EFPC's
fdXpca=pca.fd(fdX, nharm=p)
eigenvalues=fdXpca$values; scoresX=fdXpca$scores
# jth column of scoresX contains scores of the jth EFPC
harmonicsX=fdXpca$harmonics # extract the EFPC's

varianceprop=fdXpca$varprop #proportion of variance explained by the EFP's
round(varianceprop*100,0)



phi=9/4*(0:m/m)%*%t(0:m/m)
# True surface evaluated on a discrete bivariate grid
par(mar=c(1.5, 1.5, 3.5, 0.2), mfrow=c(2,2), oma = c(4, 4, 0.2, 0.2))
# 4 panels 2 rows and 2 columns arrangement
axislabelsize=1.5 # controls the size of labels of axes
axisticksize=0.8 # controls the size of ticks of axes
persp((0:m/m),(0:m/m),z=phi,cex.axis = axisticksize, cex.lab=axislabelsize, xlab="t",ylab="s", zlab=" ", theta=30,phi=30,
ticktype="detailed", main="True")

# Next we compute hat(phi)
# vivj is the matrix whose entries are products of v_j(s_k)*v_i(t_l).
# Blocks of vivj of m by m matrices represent products of # v_j(s)v_i(t) evaluated on the (m+1) by (m+1) grid
for(npc in 1:3){
  vivj=matrix(0,p*(m+1),p*(m+1))
  for(j in 1:npc){
  for(i in 1:npc){ vivj[1:(m+1)+(m+1)*(i-1),1:(m+1)+(m+1)*(j-1)] =eval.fd(evalarg=0:m/m, harmonicsX[i]) %*%t(eval.fd(evalarg=0:m/m, harmonicsX[j]))
  } }
# phip will be the estimated surface.
phip=matrix(0,m+1,m+1)
for(k in 1:(N-1)){
  temp=matrix(0,m+1,m+1)
  for(j in 1:npc){
    temp1=matrix(0,m+1,m+1)
    for(i in 1:npc){
      temp1=temp1+scoresX[k+1,i]*vivj[1:(m+1)+(m+1)*(i-1),1:(m+1)+(m+1) *(j-1)]
    }
    temp=temp+(eigenvalues[j])^(-1)*scoresX[k,j]*temp1 }
  phip=phip+temp
}

phip=(1/(N-1))*phip
if (npc==1)
  persp((0:m/m),(0:m/m), z=phip, cex.axis =axisticksize, cex.lab =axislabelsize, xlab="t",ylab="s", zlab=" ", theta=30, phi=30, ticktype="detailed", main="p=1")
else if (npc==2)
  persp((0:m/m),(0:m/m), z=phip, cex.axis = axisticksize, cex.lab=axislabelsize,
        xlab="t", ylab="s",zlab=" ",theta=30,phi=30,ticktype="detailed", main="p=2")
else if (npc==3)
  persp((0:m/m),(0:m/m), z=phip, cex.axis = axisticksize,cex.lab=axislabelsize,xlab="t", ylab="s", zlab=" ", theta=30, phi=30, ticktype="detailed", main="p=3")
}



############################################################################
# Miscellaneous
############################################################################

load("~/Desktop/Desktop/Teaching/EWHA/자료분석특론(fda)_2021_spring/FDA/slides_new/R codes/X.alcohol.RData")

X <- X.al

n=dim(X)[1]# sample size 77
ntx=dim(X)[2] # number of  time points to sample from each predictor 256
p=dim(X)[3] # number of dimension 64
nbase=15 #


# observed time points for balanced
nnt <- ntx * 5
tx <- round(seq(1, nnt, len=ntx))
tx <- matrix(tx, nrow=n, ncol=ntx, byrow=T)/ nnt


# stack X by using rbind
funcx=numeric()
for (k in 1:p){
  funcx=rbind(funcx,(X[,,k]))
}

dimen <- dim(funcx)[[2]]
size <- dim(funcx)[[1]]

# produce splines functions defined on interval [0,1]
databasis=create.fourier.basis(rangeval=c(0,1),nbasis=nbase)


# In functional data
x.fd=Data2fd(tx[1,],t(funcx),databasis)

xf <- t(x.fd$coefs)

# draw dataset
x.draw=Data2fd(tx[1,],t(X[,,1]),databasis)
plot(x.draw)
title("Alcoholic")
