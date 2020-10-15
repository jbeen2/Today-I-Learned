# 1 : AR(2) 

set.seed(1) ; n=200 ; a=rnorm(n) 
yt = function(phi){
  y = c() ; y[1] = 0 ; y[2] = 0
  for (t in 3:n) y[t] = (1+phi)*y[t-1] - phi*y[t-2] + a[t]
  return(y)
}

# 1) phi = 0 (단위근 1개)
y1 <- yt(0) 
plot(y1, type="l") # 시도표 
acf(y1) # SACF 

# 2) phi = 0.8 (단위근 1개)
y2 <- yt(0.8) 
plot(y2, type='l') # 시도표 
acf(y2) # SACF 

# 3) phi = 1 (단위근 2개)
y3 <- yt(1)
plot(y3, type='l') # 시도표 
acf(y3) # SACF 



### 2 
Edu_Ent.data <- read.table("EX3_6_Edu_Ent.txt", header=F)
y = Edu_Ent.data$V1
plot(log(y), type="o")



### 3 
library(forecast)
aic = c() 
for(p in 1:10){
  ar.fit = Arima(log(y), order=c(p,0,0))
  aic[p] = ar.fit$aic
}
plot(aic, type="o")

aic
aic.order = which.min(aic) ; aic.order


library(fUnitRoots)
adfTest(log(y), type='ct', lags=6)
