##### Brownian motion with BOA #####

library(expm); library(fda)


### PCA with Brownian Motion ###
N<-100
N_times<-500
Times<-(1:N_times)/N_times
M=200
basis<-create.bspline.basis(c(0,1),nbasis=M,norder=4)

Z<-matrix(rnorm(N*N_times),nrow=N,ncol=N_times)
BM<-t(apply(Z,1,cumsum)/sqrt(N_times))
# now expand BM functions using bsplines
BM_f<-Data2fd(Times,t(BM),basisobj=basis)
BM_pca<-pca.fd(BM_f,nharm = 6)

par(mfrow=c(1,2))
plot(BM_pca$varprop,ylab="Explained Variance per PC"); 
plot(cumsum(BM_pca$varprop),ylab="Cummulative Explained Variance")

par(mfrow=c(2,2))
for(i in 1:4){
  plot(BM_pca$harmonics[i],main=paste("PC",i,sep=""))
}


### Analysis of BOA Data ###
### Cummulative Log Returns ###
BOA<-read.table("BOA.txt",header=TRUE)
Dates<-dimnames(BOA)[[1]]
Times<-dimnames(BOA)[[2]]

N<-dim(BOA)[1]
M<-dim(BOA)[2]
T<-seq(0,6.5,length=M)

BOA<-data.matrix(BOA)
log_BOA<-log(BOA) - matrix(log(BOA)[,1],nrow=N,ncol=M)

bspline_basis<-create.bspline.basis(rangeval=c(0,6.5),norder=4,nbasis=200)
log_BOA_f<-Data2fd(T,t(log_BOA),basisobj = bspline_basis)
plot(log_BOA_f[1:10],xlab="",ylab="",lwd=1.5)

# There is an outlier #
plot(log_BOA_f)
outlier<-which(abs(log_BOA)==max(abs(log_BOA)),arr.ind=TRUE)[[1]]
log_BOA<-log_BOA[-outlier,]
N<-dim(log_BOA)[1]
M<-dim(log_BOA)[2]

# Examine the mean and variance #
Y_f<-Data2fd(T,t(log_BOA),basisobj = bspline_basis)

muhat<-mean.fd(Y_f)
sdhat<-sd.fd(Y_f)
par(mfrow=c(1,2))
plot(muhat)
plot(sdhat)

# Add some standard errors to mean #
SE_hat_U<-sdhat # create upper SE bound
SE_hat_L<-sdhat # create lower SE bound
SE_hat_U$coefs<-2*SE_hat_U$coefs/sqrt(N) + muhat$coefs
SE_hat_L$coefs<- -2*SE_hat_L$coefs/sqrt(N) + muhat$coefs
plot.fd(SE_hat_U,ylim=c(-0.002,0.002),col='red',lty=2,xlab="",ylab="")
plot.fd(SE_hat_L,add=TRUE,col='red',lty=2)
plot.fd(muhat,add=TRUE)

# Now lets examine the PCs
Y_pca<-pca.fd(Y_f,nharm = 6)
Y_pca$varprop; cumsum(Y_pca$varprop)
par(mfrow=c(2,2))
for(i in 1:4){
  plot(Y_pca$harmonics[i],xlab="",ylab="",main=paste("PC",i,sep=""))
}

plot(Y_pca$harmonics)
# Creat Brownian Motion

N=50
W.mat=matrix(0,ncol=N, nrow=10000)
times = (1:10000)/10000
for(n in 1:N){W.mat[,n]=cumsum(rnorm(10000))/100}
B.25.basis <- create.bspline.basis(rangeval=c(0,1), norder=4,nbasis=500)
W.fd <- Data2fd(times, W.mat, B.25.basis)
plot(W.fd, ylab="", xlab="", col="grey", lyt=1)

B.pca <- pca.fd(W.fd,nharm=6)

plot(B.pca$harmonics)

# or

N<-100
N_times<-500
Times<-(1:N_times)/N_times
M=200
basis<-create.bspline.basis(c(0,1),nbasis=M,norder=4)

Z<-matrix(rnorm(N*N_times),nrow=N,ncol=N_times)
BM<-t(apply(Z,1,cumsum)/sqrt(N_times))
# now expand BM functions using bsplines
BM_f<-Data2fd(Times,t(BM),basisobj=basis)
BM_pca<-pca.fd(BM_f,nharm = 6)
plot(BM_pca$harmonics)



# Theoretical PCs of Brownian Motion
t<-seq(0,1,length=100)
e1<-sqrt(2)*sin(.5*pi*t)
e2<-sqrt(2)*sin(1.5*pi*t)
e3<-sqrt(2)*sin(2.5*pi*t)
e4<-sqrt(2)*sin(3.5*pi*t)
plot(t,e1,type="l",ylim=c(-1.5,1.5),xlab="",ylab="")
points(t,e2,type="l",col=2)
points(t,e3,type="l",col=3)
points(t,e4,type="l",col=4)

# Covariance Function
Chat<-var.fd(Y_f)
eval_times<-seq(0,6.5,length=50) 
Chat_M<-eval.bifd(eval_times,eval_times,Chat) 
col <- colorRampPalette(c("blue", "purple", "red"))(n = 299)
par(mar=c(2,2,.5,.5), xaxs = "i", yaxs = "i")
image(eval_times,eval_times,Chat_M,col=col,xlab="",ylab="")

persp(eval_times,eval_times,Chat_M)