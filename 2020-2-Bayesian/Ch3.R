# Ch.3 Probability Distribution 

# 4 & 5  
p_theta = function(p) p^3 * (1-p)^2 * (1/9)
p_all = c(1:9)*0.1
p_theta(p_all)
p_theta(p_all)[5] / sum(p_theta(p_all))

post = c() 
for (i in 1:9){ 
  post[i] = p_theta(p_all)[i] / sum(p_theta(p_all))
  }
plot(post, type="o")

p_theta(p_all)[5] * 9  # conditional probability 
sum(p_theta(p_all))    # marginal probability 


# 6 
# B(10, 0.7)
hist(x <- rbinom(1000, 10, 0.7), prob=T, xlab="n=10")
curve(dnorm(x, 10*0.7, sqrt(10*0.7*0.3)), from=0, to=10, type="l", add=T)

# B(100, 0.7)
hist(x <- rbinom(1000, 100, 0.7), prob=T, xlim=range(55,85))
curve(dnorm(x, 100*0.7, sqrt(100*0.7*0.3)), from=55, to=85, type="l", add=T)
