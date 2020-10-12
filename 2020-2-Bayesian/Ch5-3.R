n=40 ; x=15 ; m=10 
a=b=1 


# true predictive dist 
z = c(0:10)
pred.z = factorial(m) / (factorial(z) * factorial(m-z)) * beta(a+x+z, b+n-x+m-z) / beta(a+x, b+n-x) 
pred.z

plot(z, pred.z, type="h") # histogram barplot


mean.z = sum(z*pred.z) 
var.z = sum(z**2 * pred.z) - mean.z**2
mean.z ; var.z


# Monte Carlo Simulation 
N = 10000
theta = rbeta(N, a+x, b+n-x)
z.sample = rbinom(N, m, theta) # 1:1 correspond 
head(theta) ; head(z.sample)

plot(table(z.sample / N), type='h')
mean(z.sample) ; var(z.sample)
