##### Lecture 11 Part 1 #####
library(fda)
library(refund)
set.seed(2016)


par(mar=c(4,4,1,1))
Y = CanadianWeather$dailyAv[,,1]
times = seq(0,1,length=365)
matplot(times,Y,type="l",ylab="Temp C",xlab="Day of Year",xlim=c(0,1),ylim=c(-40,25))
vl = cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30,31))/365
m_pts = (vl[-1] + vl[-13])/2
mts<-c("J","F","M","A","M","J","J","A","S","O","N","D")
abline(v=vl,lty=2,col="grey")
text(x=m_pts,y=-39,labels=mts,cex=1.5)


library(maps)
place <- CanadianWeather$place
coordinates <- CanadianWeather$coordinates
coordinates <- coordinates[,2:1]
coordinates[,1] <- - coordinates[,1] # longitude, latitude 위치 (x,y)

# 지도 시각화 
map('world',ylim=c(42,78),xlim = c(-170,-40))
map('world',region="Canada",add=TRUE)
points(coordinates,col='black',pch=19)
i.0 <- which(place=="Calgary")
coord.0 <- coordinates[i.0,]
points(coord.0[1],coord.0[2], pch= 15,cex=1.5,col='red')

# Ireland 지도 그리기 
library(gstat)
data(wind)
head(wind,n=3)
length(unique(wind$year))

par(mar=c(2,2,1,1))
library(sp) # char2dms
lat = as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
lon = as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))
coordinates <- cbind(lon,lat)
coordinates[,1] <- coordinates[,1]
map('world',region="Ireland")
#map('world',ylim=c(51,56),xlim = c(-10.5,-5.5))
points(coordinates,col='red',pch=19)




##### Lecture 11 Part 2 #####
library(fda)
library(refund)
set.seed(2016)

# geofd, geoR : krigging 
library(fda); library(fda.usc) 
library(geofd); library(maps)
library(geoR)
data("CanadianWeather")
dailyAv <- CanadianWeather$dailyAv
dim(dailyAv)
Temperature   <- dailyAv[,,1]
Precipitation <- dailyAv[,,2] # 강수량  
log10Precip   <- dailyAv[,,3]
place <- CanadianWeather$place
coordinates <- CanadianWeather$coordinates # 위치 : lat, lon 
#Swap lat and long
coordinates <- coordinates[,2:1]
#Flip oreientation of long 
coordinates[,1] <- -coordinates[,1]


geo.dist.35 <- dist(coordinates) # 위도와 경도를 바탕으로 각 지점마다  euclidean distance 구하기 
Day <- 1:365
n<-dim(Temperature)[2]  
nt <- dim(Temperature)[1]
i.0 <- which(place=="Calgary")
coord.0 <- coordinates[i.0,] # Calgary data 제거 -> 뺀 데이터로 Calgary 날씨 찾기 
Tempe.34 <- Temperature[,-i.0]
coord.34 <- coordinates[-i.0,]


par(mar=c(2,2,1,1))
map('world',region="canada")
points(coordinates,col='blue',pch=19)
points(coord.0[1],coord.0[2],pch=15,cex=1.5,col="red")


# 날씨는 cycle 이 있으므로, fourier basis! 
# 함수형 데이터로 변환 
K <- 99
fourier.basis <- create.fourier.basis(rangeval=
                                        range(Day),nbasis=K)
temp.fd.Fb.34 <- Data2fd(argvals=Day, y=Tempe.34, 
                         basisobj=fourier.basis)
temp.fd.Fb.0 <- Data2fd(argvals=Day, y=
                          Temperature[,i.0],basisobj=fourier.basis)


# Calgary 를 제외한 데이터를 바탕으로, Calgary 데이터 예측하기 
par(mar=c(4,4,1,1))
plot(temp.fd.Fb.34,col="grey",
     xlab="Day",ylab="Temperature (degrees C)",
     main="Average daily temperatures")
lines(temp.fd.Fb.0,lwd=2)


# computing L2 norms between functions, 
#  using Fourier basis expansions
# Note if the basis is not orthogonal, 
#  then a weight matrix is needed.
L2norm.Fb.34 <- dist(t(temp.fd.Fb.34$coefs))^2  # coef 끼리의 distance 만들기 

# Calculating the empiricial trace bin variogram
# Fitted variogram assumes an exponential covariance.
emp.trace.vari.34 <- trace.variog(coords=coord.34, 
                                  L2norm=as.matrix(L2norm.Fb.34), bin=TRUE) # semi variogram 

# fitting an exponential vriogram
sigma2.0 <- quantile(emp.trace.vari.34$v, 0.75) # p.11 C(h) = sigma^2 * phi(h)
phi.0    <- quantile(emp.trace.vari.34$Eu.d, 0.75)
fit.vari.34 <- variofit(emp.trace.vari.34, 
                        ini.cov.pars=c(sigma2.0,phi.0),cov.model="exponential") # parameter 구하는 것  


# variogram value (각 point 들의 variance ex. vancouver-London, vancouver-Sydney, ...)
par(mar=c(2,2,1,1))
plot(as.dist(emp.trace.vari.34$Eu.d),L2norm.Fb.34,col="grey",
     xlab="Geographical distances", ylab="L2 distances",
     main="Empirical variogram")
points(emp.trace.vari.34$u,emp.trace.vari.34$v,col="black",pch
       =19) # variogram 들이 비슷한 것 끼리 묶어놓은 것 (Clustering)
lines(fit.vari.34,col="black",lwd=2)  # Fitted Variogram 
legend("topleft", c("Variogram cloud",  "Binned variogram",  "
    Fitted variogram"),
       col=c(8,1,1),  lwd=c(-1,-1,2), pch=c(1,19,-1) )


# krigging \ p.21 sigma11 inverse * sigma12 -> 여기에서 sigma11 구하는 것 
hat.C.34 <- cov.spatial(emp.trace.vari.34$Eu.d,  # 모든 도시들 사이의 euclidean distance 
                        cov.model= fit.vari.34$cov.model,  
                        cov.pars=fit.vari.34$cov.pars) # fit.vari.34 parameter setting 
geo.dist.0.34 <- as.matrix(geo.dist.35)[-i.0,i.0] # Calgary 와 다른 도시들 사이의 euclidean distance 
hat.C.0 <- cov.spatial(geo.dist.0.34,
                       cov.model= fit.vari.34$cov.model,
                       cov.pars=fit.vari.34$cov.pars)

# kriging weights
w.k <- solve(hat.C.34,hat.C.0) # w11 inverse w12 
# mean est weights
inv.hat.C.34 <- solve(hat.C.34)
ones<-matrix(1,nrow=34,ncol=1)
w.m <- inv.hat.C.34%*%ones 
w.m<- w.m/c(t(ones)%*%inv.hat.C.34%*%ones) 
sum(w.m); sum(w.k)


# y : 도시들 사이의 weight / x : 도시들 간의 거리 
# krigging 에서 도시 간의 거리가 멀면 weight 작아짐 
par(mar=c(2,2,3,1))
plot(geo.dist.0.34,w.k)
abline(h=0,lty=2)


par(mar=c(2,2,0,0))
mean_coef<-temp.fd.Fb.34$coefs%*%w.m
mean_f<-fd(coef=mean_coef,basis=fourier.basis)
plot(temp.fd.Fb.34,col="grey")
plot(mean_f,lwd=2,add=TRUE) # weight 를 줘서 구한 mu hat 
# plot(mean(temp.fd.Fb.34),lwd=2,add=TRUE,lty=2,col='red') # naive mean 


# p.41 matrix multiplication 
a = 1 - sum(w.k)
tmp_coef <- temp.fd.Fb.34$coef%*%w.k # coef 에 krigging 한 것 곱해서 
tmp_fun <- fd(coef=tmp_coef,basis=fourier.basis) # 함수형 데이터 만들기 
Xhat<-a*mean_f + tmp_fun

par(mar=c(2,2,1,1))
plot(temp.fd.Fb.34,col="grey") # 43개의 grey area 
plot(temp.fd.Fb.0,add=TRUE,lwd=2) # 실제 Calgary 날씨 
plot(Xhat,add=TRUE,lty=2,lwd=2,col="red") # 우리가 생각한 Calgary 날씨 


library(gstat); library(sp)
data(wind)
lat = as.numeric(char2dms(
  as.character(wind.loc[["Latitude"]])))
lon = as.numeric(char2dms(
  as.character(wind.loc[["Longitude"]])))
coordinates <- cbind(lon,lat)  # 12개 도시, (lon, lat)
coordinates[,1] <- coordinates[,1]
coordinates<-coordinates[c(5,1,12,9,4,6,11,3,7,10,2,8),]  # reordering 
geo.dist.12 <- dist(coordinates)  # wind distance 
dim(wind) # 6574 15 : 365 일 데이터가 아닌, 18년치 데이터 
Times<-wind[,1:3]
X <- wind[,-(1:3)]


# Ireland Plot 
par(mar=c(2,2,1,1))
map('world',region="Ireland")
points(coordinates,col='blue',pch=19)
points(coordinates[7,1],coordinates[7,2],
       pch=15,cex=1.5,col="red")


years<-unique(wind$year)
K <- 99
fourier.basis <- create.fourier.basis(rangeval=
                                        c(0,1),nbasis=K)


# 다년차 데이터는 전부 다 쓰는 것 보다, 1년치씩 handling 하는 것이 좋음! 
coef_mat = matrix(0,K,12) # 도시 12개의 cov matrix 
for(i in 1:length(years)){
  year<-years[i]
  X.tmp<-as.matrix(X[wind$year==year,]) # 61년치 데이터 다 빼서 
  ndays<-dim(X.tmp)[1]
  pts<-seq(0,1,length=ndays)
  F.tmp<-Data2fd(argvals=pts,y=X.tmp,fourier.basis) # data2fd 로 만들고 
  coef_mat = coef_mat + coef(F.tmp) # 저장 
}
coef_mat<-coef_mat/18 
X.fd<-fd(coef_mat,fourier.basis)



# 7번째 도시 
temp.fd.11 <- X.fd[-7]
temp.fd.0 <- X.fd[7]
coord.11<-coordinates[-7,]
coord.0<-coordinates[7,]

# semi variable 
L2norm.11 <- dist(t(temp.fd.11$coefs))^2
emp.trace.vari.11 <- trace.variog(coords=coord.11, 
                                  L2norm=as.matrix(L2norm.11))

sigma2.0 <- quantile(emp.trace.vari.11$v, 0.75)
phi.0    <- quantile(emp.trace.vari.11$Eu.d, 0.75)
fit.vari.11 <- variofit(emp.trace.vari.11, 
                        ini.cov.pars=c(sigma2.0,phi.0),cov.model="matern")


# variogram 
par(mar=c(2,2,1,1))
plot(as.dist(emp.trace.vari.11$Eu.d),L2norm.11,col="grey",
     xlab="Geographical distances", ylab="L2 distances",
     main="Empirical variogram",xlim=c(0,5))
#points(emp.trace.vari.11$u,emp.trace.vari.11$v,col="black",pch=19)
lines(fit.vari.11,col="black",lwd=2)
legend("topleft", c("Variogram cloud","Fitted variogram"),
       col=c(8,1),  lwd=c(-1,2), pch=c(1,-1) )


i.0 <- 7
#coord.0 <- coordinates[i.0,]
#Tempe.11 <- Temperature[,-i.0]
#coord.11 <- coordinates[-i.0,]

hat.C.11 <- cov.spatial(emp.trace.vari.11$Eu.d,
                        cov.model= fit.vari.11$cov.model,
                        cov.pars=fit.vari.11$cov.pars)
geo.dist.0.11 <- as.matrix(geo.dist.12)[-i.0,i.0]
hat.C.0 <- cov.spatial(geo.dist.0.11,
                       cov.model= fit.vari.11$cov.model,
                       cov.pars=fit.vari.11$cov.pars)

# kriging weights
w.k <- solve(hat.C.11,hat.C.0)
# mean est weights
inv.hat.C.11 <- solve(hat.C.11)
ones<-matrix(1,nrow=11,ncol=1)
w.m <- inv.hat.C.11%*%ones 
w.m<- w.m/c(t(ones)%*%inv.hat.C.11%*%ones) 
sum(w.m); sum(w.k)


par(mar=c(2,2,3,1))
plot(geo.dist.0.11,w.k)
abline(h=0,lty=2)

par(mar=c(2,2,0,0))
mean_coef<-temp.fd.11$coefs%*%w.m
mean_f<-fd(coef=mean_coef,basis=fourier.basis)
plot(temp.fd.11,col="grey")
plot(mean_f,lwd=2,add=TRUE) 
#plot(mean(temp.fd.11),lwd=2,add=TRUE,lty=2,col='red') 

a = 1 - sum(w.k)
tmp_coef <- temp.fd.11$coef%*%w.k
tmp_fun <- fd(coef=tmp_coef,basis=fourier.basis)
Xhat<-a*mean_f + tmp_fun

par(mar=c(2,2,1,1))
plot(temp.fd.11,col="grey")
plot(temp.fd.0,add=TRUE,lwd=2) # 실제
plot(Xhat,add=TRUE,lty=2,lwd=2,col="red")  # 예측 

