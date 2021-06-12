library(fda)
library(refund)
library(geofd)
library(geoR)
library(maps)
library(dplyr)
library(tidyverse)
library(funHDDC)

set.seed(1234)
setwd("~/workdir/2021-1-FunctionalDataAnalysis/FinalProject")

# ======== Data Load ======== 
# 1. Temperature Data (high)
temp1991_high <- read.csv("temp1991_high_without_island.csv", header=T, encoding = "UTF-8", fileEncoding = "CP949")
temp1991_high <- temp1991_high %>% column_to_rownames(., var = "date")
temp1991_high <- as.matrix(temp1991_high)

# 2. Temperature Data (low)
temp1991_low <- read.csv("temp1991_low_without_island.csv", header=T, encoding = "UTF-8", fileEncoding = "CP949")
temp1991_low <- temp1991_low %>% column_to_rownames(., var = "date")
temp1991_low <- as.matrix(temp1991_low)

# 3. Sun : 일조시간 
sun1991 <- read.csv("sun1991.csv", header=T, encoding = "UTF-8", fileEncoding = "CP949")
sun1991 <- sun1991 %>% column_to_rownames(., var = "date")
sun1991 <- as.matrix(sun1991)

# 4. Location Data 
loc <- read.csv("location.csv", header=T, encoding = "UTF-8", fileEncoding = "CP949")
loc <- loc %>% column_to_rownames(., var = "지점명")
loc <- loc[,2:1] # swap to (lat, lon)
rownames(loc)[rownames(loc) == "대구(신암)"] <- "대구.신암."
rownames(loc)[rownames(loc) == "전주(완산)"] <- "전주.완산."
loc <- loc[colnames(temp1991_high),] 
loc <- as.matrix(loc) ; geo.dist <- dist(loc)  # euclidean distance 

# 5. Visualize  
par(mar=c(2,2,1,1))
maps::map('world', region='South Korea')  
points(loc, col='blue', pch=19)


# ======== I. Temporal Spatial FDA ========
predict_temp <- function(tempdata, i){
  # 1. Setting 
  temp1991_high <- tempdata
  n <- dim(temp1991_high)[2] ; nt <- dim(temp1991_high)[1] ; Day <- 1:366
  place <- rownames(loc)
  
  # 2. Select Location 
  i.0 <- which(place==i)  
  loc.0 <- loc[i,] 
  temp1991_high.5 <- temp1991_high[,-i.0]
  
  loc.5 <- loc[-i.0,] 

  # 3. Visualize Point 
  par(mar=c(2,2,1,1))
  maps::map('world', region='South Korea') 
  points(loc,col='blue',pch=19)
  points(loc.0[1],loc.0[2],pch=15,cex=1.5,col="red")
  legend("bottomright", c("rest",  "place chose to predict"),
         col=c("blue","red"), pch=c(16,15) )
  
  # 4. To Functional Data 
  K <- 99
  fourier.basis <- create.fourier.basis(rangeval=range(Day), nbasis=K)
  fd.temp1991_high.5 <- Data2fd(argvals=Day, y=temp1991_high.5, basisobj=fourier.basis)
  temp.fd.Fb.0 <- Data2fd(argvals=Day, y= temp1991_high.5[,i.0], basisobj=fourier.basis)

  # 5. Plot 
  par(mar=c(4,4,1,1))
  plot(fd.temp1991_high.5,col="grey", xlab="Day", ylab="Temperature (degrees C)", main="Average daily temperatures")
  lines(temp.fd.Fb.0,lwd=2)

  # 6. Variogram
  L2norm.Fb.0 <- dist(t(fd.temp1991_high.5$coefs))^2
  emp.trace.vari.34 <- trace.variog(coords=loc.5, L2norm=as.matrix(L2norm.Fb.0), bin=TRUE) # semi variogram

  sigma2.0 <- quantile(emp.trace.vari.34$v, 0.75) # p.11 C(h) = sigma^2 * phi(h)
  phi.0    <- quantile(emp.trace.vari.34$Eu.d, 0.75)
  fit.vari.34 <- variofit(emp.trace.vari.34, ini.cov.pars=c(sigma2.0,phi.0), cov.model="exponential")

  par(mar=c(2,2,1,1))
  plot(as.dist(emp.trace.vari.34$Eu.d),L2norm.Fb.0,col="grey",
       xlab="Geographical distances", ylab="L2 distances",
       main="Empirical variogram")
  points(emp.trace.vari.34$u,emp.trace.vari.34$v,col="black",pch=19)
  lines(fit.vari.34,col="black",lwd=2)  # Fitted Variogram
  legend("topleft", c("Variogram cloud", "Binned variogram", "Fitted variogram"),
         col=c(8,1,1),  lwd=c(-1,-1,2), pch=c(1,19,-1) )

  # 7. Kriging
  hat.C.34 <- cov.spatial(emp.trace.vari.34$Eu.d,
                          cov.model= fit.vari.34$cov.model,
                          cov.pars=fit.vari.34$cov.pars)
  geo.dist.0.34 <- as.matrix(geo.dist)[-i.0,i.0]
  hat.C.0 <- cov.spatial(geo.dist.0.34,
                         cov.model= fit.vari.34$cov.model,
                         cov.pars=fit.vari.34$cov.pars)

  w.k <- solve(hat.C.34,hat.C.0) # kriging weights
  inv.hat.C.34 <- solve(hat.C.34)
  ones<-matrix(1,nrow=(n-1),ncol=1)
  w.m <- inv.hat.C.34%*%ones
  w.m <- w.m/c(t(ones)%*%inv.hat.C.34%*%ones)
  sum(w.m); sum(w.k)

  par(mar=c(2,2,3,1))
  plot(geo.dist.0.34,w.k)
  abline(h=0,lty=2)

  par(mar=c(2,2,0,0))
  mean_coef<-fd.temp1991_high.5$coefs%*%w.m
  mean_f<-fd(coef=mean_coef,basis=fourier.basis)

  a = 1 - sum(w.k)
  tmp_coef <- fd.temp1991_high.5$coef%*%w.k
  tmp_fun <- fd(coef=tmp_coef,basis=fourier.basis)
  Xhat <- a*mean_f + tmp_fun

  # 8. Result
  par(mar=c(2,2,1,1))
  plot(fd.temp1991_high.5,col="grey", 
       main = paste0(i)) # Cluster
  lines(temp.fd.Fb.0,lwd=2,col="black") # Actual
  lines(mean_f,lty=2,lwd=2,col="red") # Predict
  lines(Xhat,lty=2,lwd=2,col="blue") # Predict with Kriging
  legend("topleft", c("real value",  "prediction",  "prediction with kriging"),
         col=c("black","red","blue"), lwd=c(2,2,2), lty=c(1,2,2))

}


# ======== Result ========
# [1] 예측이 잘 되었던 지역  
## 1) 관측소가 많은 지역  
predict_temp(temp1991_low, "기상청") 
predict_temp(temp1991_low, "해운대") 

## 2) 내륙 
predict_temp(temp1991_high, "남원")

## 3) kriging 시 예측력 향상 
predict_temp(temp1991_high, "군산")


# [2] 예측력이 좋지 않았던 지역 
## 1) 산맥 
predict_temp(temp1991_low, "추풍령")

## 2) 해안가 
predict_temp(temp1991_low, "포항")
predict_temp(temp1991_high, "통영")
predict_temp(temp1991_high, "울산")

## 3) kriging 시 예측력이 떨어지는 지역 
predict_temp(temp1991_high, "밀양")
predict_temp(temp1991_low, "태백")
predict_temp(temp1991_low, "태안")



# ======== II. Multivariate Clustering ========
# 1. Data Preprocessing & Filtering 
a <- intersect(colnames(temp1991_high), colnames(sun1991))
temp1991_high <- temp1991_high[,a] ; sun1991 <- sun1991[,a] 
loc <- loc[a,] ; geo.dist <- dist(loc)

# 2. To Functional Data 
K <- 99 ; Day <- 1:366
fourier.basis <- create.fourier.basis(rangeval=range(Day), nbasis=K)
fd.temphigh1991 <- Data2fd(argvals=Day, y=temp1991_high, basisobj=fourier.basis)
fd.sun1991 <- Data2fd(argvals=Day, y=sun1991, basisobj=fourier.basis)

# 3. Multivariate Clustering 
res.multi<-funHDDC(list(fd.temphigh1991, fd.sun1991), K=5, init="kmeans", threshold=0.2)
res.table<-as.data.frame(x = as.factor(res.multi$class), row.names = colnames(temp1991_high))
colnames(res.table) <- "class"

res.1 <- subset(res.table, class==1) # [1] 군산, 추풍령, 홍천   
res.2 <- subset(res.table, class==2) # [2] 대관령, 태백 (2개) 
res.3 <- subset(res.table, class==3) # [3] 고창, 대구.신암. , 합천 
res.4 <- subset(res.table, class==4) # [4] 강릉, 동해, 속초, 포항, 영덕, 군산 
res.5 <- subset(res.table, class==5) # [5] 거제, 남해, 목포, 여수, 통영  


# ======== Predict Temperature ========
predict_temp_cluster <- function(tempdata, cluster, i){
  # 1. Select Cluster 
  temp1991_high <- tempdata[,rownames(cluster)] 
  loc <- loc[rownames(cluster), ] ; geo.dist <- dist(loc)
  
  # 2. Setting 
  n <- dim(temp1991_high)[2]  ; nt <- dim(temp1991_high)[1] ; Day <- 1:366
  place <- rownames(loc)
  
  # 3. Select Location 
  i.0 <- which(place==i)  
  loc.0 <- loc[i,] 
  temp1991_high.5 <- temp1991_high[,-i.0]
  
  loc.5 <- loc[-i.0,] 
  
  # 4. Visualize Point 
  par(mar=c(2,2,1,1))
  maps::map('world', region='South Korea') 
  points(loc,col='blue',pch=19)
  points(loc.0[1],loc.0[2],pch=15,cex=1.5,col="red")
  legend("bottomright", c("rest",  "place chose to predict"),
         col=c("blue","red"), pch=c(16,15) )
  
  # 5. To Functional Data 
  K <- 99
  fourier.basis <- create.fourier.basis(rangeval=range(Day), nbasis=K)
  fd.temp1991_high.5 <- Data2fd(argvals=Day, y=temp1991_high.5, basisobj=fourier.basis)
  temp.fd.Fb.0 <- Data2fd(argvals=Day, y= temp1991_high.5[,i.0], basisobj=fourier.basis)
  
  # 6. Plot 
  par(mar=c(4,4,1,1))
  plot(fd.temp1991_high.5,col="grey", xlab="Day", ylab="Temperature (degrees C)", main="Average daily temperatures")
  lines(temp.fd.Fb.0,lwd=2)
  
  # 7. Variogram 
  L2norm.Fb.0 <- dist(t(fd.temp1991_high.5$coefs))^2 
  emp.trace.vari.34 <- trace.variog(coords=loc.5, L2norm=as.matrix(L2norm.Fb.0), bin=TRUE) # semi variogram 
  
  sigma2.0 <- quantile(emp.trace.vari.34$v, 0.75) # p.11 C(h) = sigma^2 * phi(h)
  phi.0    <- quantile(emp.trace.vari.34$Eu.d, 0.75)
  fit.vari.34 <- variofit(emp.trace.vari.34, ini.cov.pars=c(sigma2.0,phi.0), cov.model="exponential")  
  
  par(mar=c(2,2,1,1))
  plot(as.dist(emp.trace.vari.34$Eu.d),L2norm.Fb.0,col="grey",
       xlab="Geographical distances", ylab="L2 distances",
       main="Empirical variogram")
  points(emp.trace.vari.34$u,emp.trace.vari.34$v,col="black",pch=19) 
  lines(fit.vari.34,col="black",lwd=2)  # Fitted Variogram 
  legend("topleft", c("Variogram cloud", "Binned variogram", "Fitted variogram"),
         col=c(8,1,1),  lwd=c(-1,-1,2), pch=c(1,19,-1) )
  
  # 8. Kriging
  hat.C.34 <- cov.spatial(emp.trace.vari.34$Eu.d,
                          cov.model= fit.vari.34$cov.model,
                          cov.pars=fit.vari.34$cov.pars)
  geo.dist.0.34 <- as.matrix(geo.dist)[-i.0,i.0]
  hat.C.0 <- cov.spatial(geo.dist.0.34,
                         cov.model= fit.vari.34$cov.model,
                         cov.pars=fit.vari.34$cov.pars)

  w.k <- solve(hat.C.34,hat.C.0) # kriging weights
  inv.hat.C.34 <- solve(hat.C.34)
  ones<-matrix(1,nrow=(n-1),ncol=1)
  w.m <- inv.hat.C.34%*%ones
  w.m <- w.m/c(t(ones)%*%inv.hat.C.34%*%ones)
  sum(w.m); sum(w.k)

  par(mar=c(2,2,3,1))
  plot(geo.dist.0.34,w.k)
  abline(h=0,lty=2)

  par(mar=c(2,2,0,0))
  mean_coef<-fd.temp1991_high.5$coefs%*%w.m
  mean_f<-fd(coef=mean_coef,basis=fourier.basis)

  a = 1 - sum(w.k)
  tmp_coef <- fd.temp1991_high.5$coef%*%w.k
  tmp_fun <- fd(coef=tmp_coef,basis=fourier.basis)
  Xhat <- a*mean_f + tmp_fun

  # 9. Result 
  par(mar=c(2,2,1,1))
  plot(fd.temp1991_high.5,col="grey") # Cluster
  lines(temp.fd.Fb.0,lwd=2,col="black") # Actual
  lines(mean_f,lty=2,lwd=2,col="red") # Predict
  lines(Xhat,lty=2,lwd=2,col="blue") # Predict with Kriging
  legend("topleft", c("real value",  "prediction",  "prediction with kriging"),
         col=c("black","red","blue"), lwd=c(2,2,2), lty=c(1,2,2))

}


# ======== Result ========
predict_temp_cluster(temp1991_high, res.3, "대구.신암.")
predict_temp_cluster(temp1991_high, res.4, "속초")
predict_temp_cluster(temp1991_low, res.5, "여수")


