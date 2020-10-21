### Chapter 2. Fundamental Concepts 

# AR(1) : yt = 1+0.8yt-1+at
n = 1000 ; a = rnorm(n) ; y = rep(n,0)
y[1] = 5 
for (t in 2:n){
  y[t] = 1+0.8*y[t-1]+a[t]
}
plot(y, type="l")
mean(y) ; var(y)
acf(y) 
pacf(y)

# I(1), Random Walk : yt = yt-1 + at 
n = 1000 ; a=rnorm(n) ; y = rep(n,0)
y[1] = 0 
for (t in 2:n){
  y[t] = y[t-1] + a[t]
}
plot(y, type="l")
mean(y) ; var(y) ; mean(y[1:200]) ; mean(y[201:401])
acf(y) 
pacf(y)


# PACF 
VKOSPI.data = read.csv("VKOSPI_2004_2019.csv", header=T)
y = VKOSPI.data$VKOSPI ; n=length(y) 

y0 = c() ; y1=c() ; y2 = c() 
y0 = y[3:n]; y1 = y[2:(n-1)] ; y2=y[1:(n-2)]

var(y0) ; cov(y0, y1) ; cov(y0,y1) / var(y0) ; acf(y)[1]
lm1 = lm(y0~y1) ; lm2 = lm(y2~y1)
e1 = lm1$residuals ; e2 = lm2$residuals
cor(e1, e2)
pacf(y)[2]


# k=3 
y0 = y[5:n]; y1 = y[4:(n-1)] ; y2=y[3:(n-2)] ; y3=y[2:(n-3)] ; y4 = y[1:(n-4)]
k0 = 4

# Durbin Levinson Recursion 
rho = acf(y)$acf[2:25]
phi = matrix(0, ncol=k0, nrow=k0)
phi[1,1] = rho[1]

for (k in 1:(k0-1)){
  phi[k+1, k+1] = (rho[k+1]-sum([phi[k,1:k]*rho[k:1]]))/(1-sum(phi[k,1:k]*rho[1:k])) ; 
  for (j in 1:k) { 
    phi[k+1,j] = phi[k,j] - phi[k+1,k+1]*phi[k,k+1-j]
    }  
}
 

# KOSPI log-return 
KOSPI.data = read.csv("KOSPI_2004_2019.csv", header=T)
return = diff(KOSPI.data$KOSPI)

acf(return, 25, main="ACF")
pacf(return, 25, main="PACF")


# HAC se 
library(sandwich)
lm.return = lm(return~1)
HAC.se.mu = sqrt(vcovHAC(lm.return)) ; HAC.se.mu
sd(return) / sqrt(length(return))

lm.VKOSPI = lm(VKOSPI.data$VKOSPI~1)
HAC.se.VKOSPI = sqrt(vcovHAC(lm.VKOSPI)) ; HAC.se.VKOSPI 
sd(VKOSPI.data$VKOSPI) / sqrt(length(VKOSPI.data$VKOSPI))


# SACF 
# AR(1) : yt = phi*yt-1 + at 
set.seed(1) ; n=100 ; phi=0.8 ; a=rnorm(n) 
y = c() ; y[1] = 0 
for (t in 2:n) {
  y[t] = phi*y[t-1] + a[t]
}
acf(y)
acf(y, type=c("covariance"))  # autocovariance 
acf(y, type=c("covariance"))$acf[1]
var(y)

set.seed(1) ; n=100 ; M=200 ; rho.1 = c() ; phi=0.8
for (i in 1:M) {
  a = rnorm(n) ; y=c() ; y[1] = 0 
  for (t in 2:n) {
    y[t] = phi*y[t-1] + a[t]
  }
  rho.1[i] = acf(y, plot=FALSE)$acf[2]
}
plot(rho.1, type="l")
mean(rho.1) ; sd(rho.1)
