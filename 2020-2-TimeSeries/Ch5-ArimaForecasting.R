### Chapter 5. ARIMA Forecasting 

# AR(1) vs AR(2) vs MA(1) 
Ex.data3_2 = read.table("Criminal.txt", header=F)
y = Ex.data3_2$V1

library(forecast)

# AR(1) 
AR1.fit = Arima(y, order=c(1,0,0)) ; AR1.fit

AR1.hat = forecast(AR1.fit, h=20) ; AR1.hat
plot(AR1.hat)


# AR(2) 
AR2.fit = Arima(y, order=c(2,0,0)) ; AR2.fit

AR2.hat = forecast(AR2.fit, h=20)
plot(AR2.hat)


# MA(1) 
MA1.fit = Arima(y, order=c(0,0,1)) ; MA1.fit
MA1.hat = forecast(MA1.fit, h=20)
plot(MA1.hat)

a.hat = MA1.fit$residuals 
plot(a.hat, type="b")



# Out-of-sample 1-step ahead forecast 
n = length(y) ; m=round(0.15*n) 

# AR(1) 
e.ar1 = c() 
for (k in 1:m) {
  N = n-k ; ar1 = Arima(y[1:N], order=c(1,0,0)) ; 
  y.hat.ar1 = forecast(ar1, h=1)$mean[1] 
  e.ar1[k] = y[N+1] - y.hat.ar1
}

# AR(2)
e.ar2 = c() 
for (k in 1:m) {
  N = n-k ; ar2 = Arima(y[1:N], order=c(2,0,0)) ; 
  y.hat.ar2 = forecast(ar2, h=1)$mean[1] 
  e.ar2[k] = y[N+1] - y.hat.ar2
}

# MA(1) 
e.ma1 = c() 
for (k in 1:m) {
  N = n-k ; ma1 = Arima(y[1:N], order=c(0,0,1)) ; 
  y.hat.ma1 = forecast(ma1, h=1)$mean[1] 
  e.ma1[k] = y[N+1] - y.hat.ma1
}


# RMSE, MAE, MAPE 
sd(e.ar1) ; sd(e.ar2) ; sd(e.ma1)  # RMSE 
mean(abs(e.ar1)) ; mean(abs(e.ar2)) ; mean(abs(e.ma1)) ;  # MAE 

y.test = c()  # MAPE 
for (k in 1:m) { 
  N = n-k ; y.test[k] = y[N+1]
}
100*mean(abs(e.ar1)/y.test) ; 100*mean(abs(e.ar2)/y.test) ; 100*mean(abs(e.ma1)/y.test)



# Eventual Forecast Functions 
Beer.data = read.table("Beer.txt", header=F)
y = Beer.data$V1
plot(y, type="l")

# 단위근 0개 : ARIMA(1,0,0)
I0.fit = Arima(y, order=c(1,0,0))
plot(forecast(I0.fit, h=1000))
I0.fit

# 단위근 1개 : ARIMA(1,1,0)
I1.fit = Arima(y, order=c(1,1,0))
plot(forecast(I1.fit, h=10))
I1.fit

# 단위근 2개 : ARIMA(1,2,0)
I2.fit = Arima(y, order=c(1,2,0))
plot(forecast(I2.fit, h=10))
I2.fit
