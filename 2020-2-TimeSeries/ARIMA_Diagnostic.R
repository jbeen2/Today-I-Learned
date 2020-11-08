GDP.data=read.csv("Kor_GDP_1960_2018.csv", head=T); attach(GDP.data); y=GDP

plot(y, type="b")      # 원계열 시도표
plot(log(y), type="b") # 로그변환 시도표 

# ARIMA
library(forecast)
aic=c(); log.y = log(y)
for(p in 1:10){ar.fit = Arima(log.y, order = c(p,0,0))
aic[p] = ar.fit$aic};  which.min(aic)

# Unit Root Test 
library(fUnitRoots)
adfTest(log.y, type='ct', lags=3)

# AR fitting 
AR.4 = Arima(log.y, order=c(4,0,0)) ; AR.4

# Diagonastic Checking 
Box.test(AR.4$residuals, lag=3, type="Ljung-Box")



# diff 
d.log.y = diff(log.y)
acf(d.log.y) 
pacf(d.log.y)

# ARIMA fitting with BIC 
bic=c(); 
for(p in 1:10){ar.fit2 = Arima(d.log.y, order = c(p-1,1,0))
bic[p] = ar.fit2$bic};  which.min(bic)

arima210.fit = Arima(d.log.y, order=c(2,1,0)) ; arima210.fit

# Portmanteau 
L = c() ; Q = c() ; p = c() 
for (i in 1:4) { L[i] = 3*i ; 
Q[i] = Box.test(arima210.fit$residuals, lag=L[i] , type="Ljung-Box")$statistic 
df = L[i]-3 
p[i] = 1-pchisq(Q[i], df) }
Portmanteau = cbind(L, Q, p) ; Portmanteau

# residual's ACF, PACF 
acf(arima210.fit$residuals)
pacf(arima210.fit$residuals)


# Forecasts 
arima210.hat = forecast(arima210.fit, h=10) ; arima210.hat
plot(arima210.hat)
