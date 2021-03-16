##### Lecture 1. Part 1 #####

library(fda)
library(ggplot2)

# Berkeley Growth Study : 함수를 smooth 하게 만들어서, shape 를 본다 
matplot(growth$age,growth$hgtf,xlab="Age",ylab="Height (cm)",pch=1,type="b",cex.lab=1.25,cex.axis=1.25)


# periodic data (시작점과 끝점이 같을 때 = 주기성이 있는 데이터)
par(mar=c(4,4,1,1))
fourier_five <- create.fourier.basis(c(0,1), nbasis=5)
plot(fourier_five,lwd=2)


# splines : 자료를 쪼개서 각각의 polynoimal 로 근사시킨다 
# spline of order n : piecewise polynomial w/ degree n-1
myfun<-function(x){sin(x)}
obs<-seq(0,2*pi,length=100)
X<-myfun(obs)
knts<-c(2,3,4)
par(mfrow=c(1,3))
for(i in knts){
  breaks<-seq(0,2*pi,length=i)
  bspline_basis<-create.bspline.basis(c(0,2*pi),breaks=breaks) 
  X.f<-Data2fd(obs,X,bspline_basis)
  curve(myfun,from=0,to=2*pi,ylab="",xlab="",cex.axis=1.25, main=paste("Knots = ",i))
  plot(X.f,add=TRUE,lty=2)
  abline(v=breaks,lty=2,col='red')
}


# Bsplines (Basis Splines)
knts<-c(2,3,4)
par(mfrow=c(1,3))
for(i in knts){
  breaks<-seq(0,2*pi,length=i)
  bspline_basis<-create.bspline.basis(c(0,2*pi),breaks=breaks)
  plot(bspline_basis, main=paste("Knots = ",i))
} 


par(mar=c(4,4,1,1))
times = growth$age; GHeight = growth$hgtf
my_basis<-create.bspline.basis(c(1,18),nbasis=10)
GHeight.F<-Data2fd(times,GHeight,my_basis)  # GHeight.F : 변환된 데이터 
plot(GHeight.F)

class(GHeight.F)   # class : 'fd'
names(GHeight.F)   # names : 'coefs', 'basis', 'fdnames'


# 1. mean function 
par(mar=c(4,4,1,1))
par(mfrow=c(1,2))
plot(GHeight.F, col='grey')
plot(mean.fd(GHeight.F), lwd=4, add=TRUE)   # coef mean, 추세 선 파악 


# 2. covariance function 
library(plot3D)
par(mar=c(4,4,1,1))
GHeight_var<-var.fd(GHeight.F)  
pts<-seq(from=1,to=18,length=50)
GHeight_mat = eval.bifd(pts,pts,GHeight_var)
persp3D(pts,pts,GHeight_mat)   # cov(s,t) 크다 = Xn(s), Xn(t) 가 동시에 평균보다 훨씬 높거나/낮다  


# 3. Functional Principal Components 
par(mar=c(4,4,1,1))
GHeight_pca<-pca.fd(GHeight.F, nharm=2)
plot(GHeight_pca$harmonics)

First_PC<-round(GHeight_pca$varprop[1],digits=3)*100
Second_PC<-round(GHeight_pca$varprop[2],digits=3)*100
First_PC ; Second_PC   # first FPCs explain 89.3% , second FPCs explain 6.4% of variability 



##### Lecture 1. Part 2 #####
### Chapter 1. First Steps in the analysis of functional data 

# *** Xn(t) = sigma Cnm Bm(t) *** 

# Mean Function : Berkeley Growth 
par(mar=c(4,4,1,1))
times = growth$age; GHeight = growth$hgtf
my_basis<-create.bspline.basis(c(1,18),nbasis=10)
GHeight.F<-Data2fd(times,GHeight,my_basis)
mu.F<-mean.fd(GHeight.F)
plot(GHeight.F, col='gray') ; par(new=TRUE) ; plot(mu.F)

# Mean Function : Alt Calc 
par(mar=c(4,4,1,1),mfrow=c(1,2))
GHeight_coef<-GHeight.F$coefs
mean_coef<-rowMeans(GHeight_coef)
mu.F2<-fd(coef=mean_coef,basisobj=my_basis)
plot(mu.F,ylab="",xlab="");plot(mu.F2,ylab="",xlab="")


# Covariance Function - bifd (bivariate functional object)
GHeight_var<-var.fd(GHeight.F)
class(GHeight_var)   # "bifd" 
names(GHeight_var)   # "coefs", "sbasis", "tbasis", "bifdnames"
dim(GHeight_var$coefs)   # 10 10 


# Covariance Function - Manual 
coef_mat<-coef(GHeight.F); N<-dim(coef_mat)[2]   # N=54 
coef_mean<-rowMeans(coef_mat)   # coef 각각에 대한 mean 값 나옴 (10개)
coef_center<-sweep(coef_mat,1,coef_mean)   # sweep : 통계량 적용하는 함수 
C_hat_coef<-(coef_center%*%t(coef_center))/(N-1)
C_hat<-bifd(coef=t(C_hat_coef),sbasisobj = my_basis, 
            tbasisobj = my_basis)


par(mar=c(1,1,1,2)) #,mfrow=c(1,2))
pts<-seq(from=1,to=18,length=100)
GHeight_mat = eval.bifd(pts,pts,GHeight_var)
C_hat_mat = eval.bifd(pts,pts,C_hat)


# cov(s,t) high = Xn(t), Xn(s) 가 동시에 평균에서 많이 벗어났다 
persp3D(pts,pts,GHeight_mat,axes=FALSE,colkey=FALSE,bty="g",shade=.4) # main="Auto" 
persp3D(pts,pts,C_hat_mat,main="Manual")   # 직접 식 구현한 것 


# Functional Principal Components Analysis (FPCA)
# direction of the greatest variation, covariance 의 축을 알고 싶다 
GHeight_pc<-pca.fd(GHeight.F,nharm=3)
class(GHeight_pc)   # 'pca.fd' 
names(GHeight_pc)[1:3]; names(GHeight_pc)[4:5]

# ========= FPCA in R =========
# nharm : the number of desired pcs 
# harmonics : principal components 
# values : eigenvalues 
# scores : the coefficients in the basis expansion 
# varprop : explained variances for each pc (ex. 89.3% / 6% / 0.2% / ...)
# meanfd : mean funciton of data 
# =============================


# FPCA in R : manual 
library(expm); 
W2inv<-inprod(my_basis,my_basis);   # inner product 
Sig_c<-cov(t(coef_mat)); A<-Sig_c%*%W2inv   # cov 구해서 PCA 
e_A<-eigen(A)
names(e_A)
e_A$values[1:3]
GHeight_pc$values[1:3]

par(mar=c(1,1,1,2),mfrow=c(1,2))
v_coef<-e_A$vectors[,1:3]
norm_v<-diag(t(v_coef)%*%W2inv%*%v_coef)
v_coef<-v_coef%*%diag(1/sqrt(norm_v))
e_fun<-fd(v_coef,my_basis)
plot(e_fun); plot(GHeight_pc$harmonics)


# FPCA - scores 
AllHeight = cbind(growth$hgtm,growth$hgtf)
N1<-dim(growth$hgtm)[2];N2<-dim(growth$hgtf)[2]
Gender<-c(rep(1,times=N1),rep(2,times=N2))
AllHeight.F<-Data2fd(times,AllHeight)
AllHeight.pc<-pca.fd(AllHeight.F,nharm=2)
Scores<-AllHeight.pc$scores
plot(Scores,col=Gender,xlab="PC1",ylab="PC2",main="Plot of scores")
legend(x="topleft",col=1:2,legend=c("Male","Female"),pch=1)


# FPCA - explained Variance 
# the first p eigenfunctions explain the following proportion of variance 
par(mar=c(3,3,1,1),mfrow=c(1,2))
AllHeight.pc<-pca.fd(AllHeight.F,nharm=10)
lamb<-AllHeight.pc$values
plot(lamb[1:10],xlab="Eigenvalue"); plot(cumsum(lamb[1:10])/sum(lamb),xlab="Explained Variance")


# Diffusion Tensor Imaging
library(refund)   # DTI dataset 
par(mar=c(3,3,1,1))
Corp<-DTI$cca
drop<-unique(which(is.na(Corp),arr.ind=TRUE)[,1])
Corp<-Corp[-drop,] # Missing value
pts<-seq(0,1,length=93)
my_basis<-create.bspline.basis(c(0,1),nbasis=20)
Corp.F<-Data2fd(pts,t(Corp),my_basis)
plot(Corp.F)


# DTI mean 
par(mar=c(3,3,1,1))
plot(Corp.F,col="gray")
Corp.F.mean<-mean(Corp.F) ; par(new=TRUE)
plot(Corp.F.mean, lwd=2)

# DTI covariance
par(mar=c(3,3,1,1))
Cov_DTI<-var.fd(Corp.F)
pts<-seq(0,1,length=200)
Cov_DTI_mat<-eval.bifd(pts,pts,Cov_DTI)
persp3D(pts,pts,Cov_DTI_mat,axes=FALSE,colkey=FALSE,bty="g",col="white",shade=.5)

DTI_pc<-pca.fd(Corp.F,nharm=2)
DTI_pc$varprop   

par(mar=c(3,3,1,1))
plot(DTI_pc$harmonics)   # weighted average 

par(mar=c(3,3,1,1))
DTI_e_fun<-DTI_pc$harmonics
DTI_scores<-DTI_pc$scores
dim(DTI_scores)
approx_fun<-Corp.F.mean+DTI_e_fun[1]*DTI_scores[1,1]+ 
  DTI_e_fun[2]*DTI_scores[1,2] 
plot(Corp.F[1]); par(new=TRUE); plot(approx_fun, lty=2)



