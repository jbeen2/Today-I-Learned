### Chapter 3. Stationary Time Series Models 

# AR(1) : Zt = 1 + 0.9Zt-1 + at 
set.seed(1) ; n=250 ; phi=0.9 ; mu=10 ; a=rnorm(n) ; z=rep(n,0)
z[1] = 10 
for (t in 2:n){
  z[t] = mu*(1-phi) + phi*z[t-1] + a[t]
}
plot(z, type="l")
mean(z) ; var(z)

# Population ACF 
rho=c() ; M=25 ; rho[1] = 1 
for (k in 2:M) {
  rho[k] = phi*rho[k-1]
}
plot(rho)
acf(z)
pacf(z)


# AR(1) : Zt = 10(1-(-.65)) + (-.65)*Zt-1 + at 
set.seed(1) ; n=250 ; phi=-0.65 ; mu=10 ; a=rnorm(n) ; z = rep(n,0)
z[1] = 10 
for (t in 2:n){
  z[t] = mu*(1-phi) + phi*z[t-1] + a[t]
}
plot(z, type="l")
mean(z) ; var(z)

# PACF 
rho=c() ; M=25 ; phi=-0.65 ; rho[1] = 1 
for (k in 2:M) {
  rho[k] = phi*rho[k-1]
}
plot(rho, type="l")
acf(z)
pacf(z)



# AR(2) : (1+0.5B-0.3B^2)Zt = at 
set.seed(1) ; n=250 ; phi.1 = -.5 ; phi.2 = .3 
a = rnorm(n) ; z = c() ; z[1] = 0 ; z[2] = 0 
for (t in 3:n) {
  z[t] = phi.1*z[t-1] + phi.2*z[t-2] + a[t]
}
plot(z, type="l")
mean(z) ; var(z)

# Population ACF  
rho = c() ; M = 25 
rho[1] = 1 ; rho[2] = phi.1 / (1-phi.2)
for (k in 3:M) {
  rho[k] = phi.1 * rho[k-1] + phi.2 * rho[k-2]
}
plot(rho, type="l")
acf(z)
pacf(z)



# AR(2) : (1-0.9B+0.9B^2)Zt = at 
set.seed(1) ; n=250 ; phi.1 = .9 ; phi.2 = -.9 
a = rnorm(n) ; z = c() ; z[1] = 0 ; z[2] = 0 
for (t in 3:n) {
  z[t] = phi.1*z[t-1] + phi.2*z[t-2] + a[t]
}
plot(z, type="l")
mean(z) ; var(z)

# Population ACF  
rho = c() ; M = 25 
rho[1] = 1 ; rho[2] = phi.1 / (1-phi.2)
for (k in 3:M) {
  rho[k] = phi.1 * rho[k-1] + phi.2 * rho[k-2]
}
plot(rho, type="l")
acf(z)
pacf(z)



# MA(1) : Zt = (1-0.5B)at 
n=250 ; set.seed(1) ; a=rnorm(n) 
z=c() ; theta=0.5 ; z[1] = 0 
for (t in 2:n) { 
  z[t] = a[t] - theta*a[t-1]
  }
plot(z, type="l")
acf(z)
