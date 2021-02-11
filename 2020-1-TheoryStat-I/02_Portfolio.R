set.seed(1234)

mu <- 0.05
sd <- 0.3

z <- matrix(rnorm(1000*20, mu, sd), nrow=1000, ncol=20)

n <- 5
v0 <- 1
rf <- 0.02
rstar <- 0.04
a <- 0.5

v1<-cbind(matrix(1,nrow=1000,ncol=1),matrix(0,nrow=1000,ncol=20))
v2<-cbind(matrix(1,nrow=1000,ncol=1),matrix(0,nrow=1000,ncol=20))

for (i in 1:1000) {
  for (j in 2:21) {
    v1[i,j]<-v1[i,j-1]*(1+(1-a)*rf+a*z[i,j-1])
    v2[i,j]<-v2[i,1]*((1+rstar)^(j-1))
    p<-colMeans(ifelse(v1>v2,1,0))
  }
}

p[2:21]


# b
plot(v1[1,1:11],type="l",ylab="Vt(a)",col=1,ylim=c(0,3))
lines(v1[2,1:11],type="l",col=2)
lines(v1[3,1:11],type="l",col=3)
lines(v1[4,1:11],type="l",col=4)
lines(v1[5,1:11],type="l",col=5)


# c 
a2<-matrix(seq(0,1,by=0.1),ncol=1)
v1_2<-cbind(matrix(1,nrow=1000,ncol=1),matrix(0,nrow=1000,ncol=20))
v2_2<-cbind(matrix(1,nrow=1000,ncol=1),matrix(0,nrow=1000,ncol=20))
p2<-matrix(rep(0,11*21),nrow=11)

for (k in 1:11) {
  for (i in 1:1000) {
    for (j in 2:21) {
      v1_2[i,j]<-v1_2[i,j-1]*(1+(1-a2[k,1])*rf+a2[k,1]*z[i,j-1])
      v2_2[i,j]<-v2_2[i,1]*((1+rstar)^(j-1))
    }
  }
  p2[k,c(1:21)]<-colMeans(ifelse(v1_2>v2_2,1,0))
}

plot(a2[,1],p2[,n+1],type="b",xlab="a",ylab="Pr(Vn(a)>Vn*)")


# d 
astar<-a2[which.max(p2[,n+1]),1]
astar

RSFC<-max(p2[,n+1])
RSFC

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


simulation(0.05,0.3,0.03,5)
simulation(0.05,0.3,0.03,10)
simulation(0.05,0.3,0.04,5)
simulation(0.05,0.3,0.04,10)
simulation(0.05,0.3,0.06,5)
simulation(0.05,0.3,0.06,10)
