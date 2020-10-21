### Chapter 6. Model Identification 

# AR(p) vs MA(q) 
Birth.Death.data = read.csv("Kor_Birth_Death.csv", header=T)
attach(Birth.Death.data)

plot(Birth, type="l")
plot(Death, type="l")

d.Birth = diff(Birth) ; d.Death = diff(Death)
plot(d.Birth, type="l")
plot(d.Death, type="l")

acf(d.Birth)
pacf(d.Birth)

acf(d.Death)
pacf(d.Death)


# model identification with AIC & BIC 
library(forecast)
d.Birth.ar.aic = c() ; d.Birth.ar.bic = c() ; 

for (p in 1:10) { 
  d.Birth.ar.fit = Arima(d.Birth, order=c(p,0,0))
  d.Birth.ar.aic[p+1] = d.Birth.ar.fit$aic 
  d.Birth.ar.bic[p+1] = d.Birth.ar.fit$bic 
}

plot(d.Birth.ar.aic, type="b")
plot(d.Birth.ar.bic, type="b")


d.Death.ar.aic = c() ; d.Death.ar.bic = c() ; 

for (p in 1:10) { 
  d.Death.ar.fit = Arima(d.Death, order=c(p,0,0))
  d.Death.ar.aic[p+1] = d.Death.ar.fit$aic 
  d.Death.ar.bic[p+1] = d.Death.ar.fit$bic 
}

plot(d.Death.ar.aic, type="b")
plot(d.Death.ar.bic, type="b")
