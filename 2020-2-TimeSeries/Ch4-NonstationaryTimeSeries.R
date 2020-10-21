### Chapter 4. Nonstationary Time Series Models 

# Random Walk : (1-B)Zt = at 
set.seed(1) ; n=100 ; a=rnorm(n) ; 
z=c() ; z[1] = 0 
for (t in 2:n){
  z[t] = z[t-1]+a[t]
}
plot(z, type="l")
acf(z)

# Random Walk with drift : (1-B)Zt = theta0 + at, t*theta0 : deterministic trend 
set.seed(1) ; n=100 ; a=rnorm(n) 
z=c() ; z[1] = 0 ; theta0 = .1 
for (t in 2:n) { 
  z[t] = theta0 + z[t-1] + a[t] # Z0 + t*theta0 + (at + at-1 + ... + a1)
}
plot(z, type="l")



# Exponential Smoothing : IMA(1,1) 
hotel.data = read.table("Hotel_month_2000_2005.txt")
hotel.data = data.frame(hotel = t(hotel.data))

library(forecast)
ses.fit = ses(hotel.data$hotel, alpha=0.2)
plot(ses.fit)



# Unit Root Test 
Ex.data3_4 = read.table("EX3_4_Beer.txt", header=F)
y = Ex.data3_4$V1 
plot(y, type="l")  # H0 : Stochastic Trend  vs  H1 : Deterministic Trend 
acf(y, main="ACF")
pacf(y, main="PACF")

aic = c() 
for (p in 1:10) { 
  ar.fit = Arima(y, order=c(p,0,0))
  aic[p] = ar.fit$aic
}
plot(aic, type="o") # lags=5 

library(fUnitRoots)
adfTest(y, type="ct", lags=5)