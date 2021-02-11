setwd("C:\\Users\\LG\\Desktop\\stat\\theorystats\\0331_assignment2")

# a 
library(tidyverse)

kospi<-read.csv("kospi.csv", sep=",")
kospi$lnPt <- log(kospi$Pt)

ggplot(kospi,aes(date,Pt,group=1))+geom_line()
ggplot(kospi,aes(date,lnPt,group=1))+geom_line()


# b
for (m in 2:nrow(kospi)) {
  kospi<-kospi %>% mutate(dlnPt=lnPt-lag(lnPt),rt=Pt/lag(Pt)-1)
}
kospi[1,c(4,5)]<-0

ggplot(kospi,aes(date,dlnPt,group=1))+geom_line()
ggplot(kospi,aes(dlnPt))+geom_histogram()


# c 
qqnorm(kospi$dlnPt) ; qqline(kospi$dlnPt, col=2)


# d 
muhat <- mean(kospi$dlnPt) / (1/12)
muhat

sigmahat <- sqrt(var(kospi$dlnPt) / (1/12))
sigmahat


# e 
simulation<-function(mu,sd,rstar,period) {
  
  set.seed(1234)
  z<-matrix(rnorm(1000*20,mu,sd),nrow=1000,ncol=20)
  a_new<-matrix(seq(0,1,by=0.1),ncol=1)
  v1_new<-cbind(matrix(1,nrow=1000,ncol=1),matrix(0,nrow=1000,ncol=20))
  v2_new<-cbind(matrix(1,nrow=1000,ncol=1),matrix(0,nrow=1000,ncol=20))
  p_new<-matrix(rep(0,11*21),nrow=11)
  
  for (k in 1:11) {
    for (i in 1:1000) {
      for (j in 2:21) {
        v1_new[i,j]<-v1_new[i,j-1]*(1+(1-a_new[k,1])*rf+a_new[k,1]*z[i,j-1])
        v2_new[i,j]<-v2_new[i,1]*((1+rstar)^(j-1))
      }
    }
    p_new[k,c(1:21)]<-colMeans(ifelse(v1_new>v2_new,1,0))
  }
  astar<-a_new[which.max(p_new[,period+1]),1]
  RSFC<-max(p_new[,period+1])
  return(list(astar,RSFC))
}

v0 <- 1
rf <- 0.02
rstar <- 0.04
a <- 0.5

v1<-cbind(matrix(1,nrow=1000,ncol=1),matrix(0,nrow=1000,ncol=20))
v2<-cbind(matrix(1,nrow=1000,ncol=1),matrix(0,nrow=1000,ncol=20))


simulation(muhat,sigmahat,0.04,5)
simulation(muhat,sigmahat,0.04,10)
simulation(muhat,sigmahat,0.04,20)
simulation(muhat,sigmahat,0.05,5)
simulation(muhat,sigmahat,0.05,10)
simulation(muhat,sigmahat,0.05,20)

library(knitr)

