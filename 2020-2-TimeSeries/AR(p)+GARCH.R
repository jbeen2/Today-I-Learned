# Data Load 
rv.data=read.csv("Won_Dollar_P_RV_r_2006_2012.csv", head=T)
attach(rv.data)
plot(RV, type="l")


# 1. AR(p)
library(forecast)

bic=c()
for(p in 1:10){ar.fit = Arima(RV, order = c(p,0,0))
bic[p] = ar.fit$bic};  which.min(bic)
plot(bic, type="b")

# AR fitting 
AR.6 = Arima(RV, order=c(6,0,0)) ; AR.6


# 2. Volatility Clustering
at = AR.6$residuals ; plot(at)


# 3. at^2 = autocorrelated 
acf(at^2)


# 4. AR(p) + GARCH(1,1)
library(rugarch)

spec.garch.t = ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1)),
               mean.model=list(armaOrder=c(0,0), include.mean=F), distribution.model="std")

garch.fit = ugarchfit(data = at, spec=spec.garch.t) ; garch.fit


# 5. sigma(T+1) 5% VaR 
ugarchforecast(garch.fit, n.ahead = 1)@forecast$sigmaFor


# 6. EGARCH(1,1) + normal 
spec.egarch = ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),
                         mean.model=list(armaOrder=c(0,0), include.mean=F), distribution.model="norm")

egarch.fit = ugarchfit(data = at, spec=spec.egarch) ; egarch.fit



# macOS error.. :( solution 
# install.packages("mclust", repo = 'https://mac.R-project.org')
