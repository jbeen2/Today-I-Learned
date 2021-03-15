# 2.2 Metropolis-Hastings Algorithm 

# random walk 
# 1. select next sample : left or right 로 이동하며 후보 값 선택 
theta = c(1:5) ; prob = c(1,3,8,5,3)/20 

prob.ratio = function(theta1, theta2, theta, prob) { 
  ind1 = which(theta == theta1)
  ind2 = which(theta == theta2)
  return (prob[ind1] / prob[ind2])
}

N = 50000
theta.curr = 2   # current sample 
theta.Samples = c(1:N)*0
theta.Samples[1] = theta.curr   # current sample 에서 시작 

# ------ Start Simulation ------ # 
for (iter in 1:(N-1)){
  # 0.5 보다 작으면 left sample, 0.5보다 크면 right sample 
  theta.prop = ifelse(runif(1) < 0.5, theta.curr+1, theta.curr-1)
  if (theta.prop < 1 || theta.prop > 5){
    theta.prop = theta.curr
    theta.prop = round(theta.prop, 0)}   # numeric error 방지 
  
  # p(theta.prop) / p(theta.curr) 확률로 후보값으로 이동, 
  # {1 - p(theta.prop) / p(theta.curr)} 확률로 현재값에 머무름 
  alpha.star = prob.ratio(theta.prop, theta.curr, theta, prob)
  
  # p(theta.prop) > p(theta.curr) 라면 alpha = 1 이므로 무조건 theta.prop 으로 이동 
  # p(theta.prop) < p(theta.curr) 라면 alpha < 1 이므로 theta.prop or 제자리 
  alpha = min(1, alpha.star)
  
  # unif~(0,1) 에 대해, u < alpha 이면 theta.next = theta.prop
  # u > alpha 이면 theta.next = theta.curr 
  theta.next = ifelse(runif(1) < alpha, theta.prop, theta.curr)
  theta.Samples[iter+1] = theta.next
  theta.curr = theta.next 
}
# ------ End Simulation ------ # 

hist(theta.Samples[500:50000], breaks=seq(0,5,by=1))



# 2. select next sample : 같은 확률로 sample 내의 후보 값 선택 
theta.Samples = c(1:N)*0
theta.Samples[1] = theta.curr

# ------ Start Simulation ------ # 
for (iter in 1:(N-1)){
  theta.prop = sample(theta, 1)  # 5개 theta sample 중 하나 임의로 선택  
  alpha.star = prob.ratio(theta.prop, theta.curr, theta, prob)
  alpha = min(1, alpha.star)
  theta.next = ifelse(runif(1) < alpha, theta.prop, theta.curr)
  theta.Samples[iter+1] = theta.next
  theta.curr = theta.next 
}
# ------ End Simulation ------ # 

hist(theta.Samples[500:50000], breaks=seq(0,5,by=1))



# 3-1. select next sample : 각각 다른 확률로 sample 내의 후보 값 선택 
theta.Samples = c(1:N)*0
theta.Samples[1] = theta.curr

# ------ Start Simulation ------ # 
for (iter in 1:(N-1)){
  # 5개 theta 중 하나 각각 다른 확률로 선택 
  prob.prop = c(.1, .1, .2, .3, .3)   # (1,2,3,4,5) 선택 될 확률 
  theta.prop = sample(theta, 1, prob=prob.prop)     
  
  alpha.star = prob.ratio(theta.prop, theta.curr, theta, prob)
  alpha = min(1, alpha.star)
  theta.next = ifelse(runif(1) < alpha, theta.prop, theta.curr)
  theta.Samples[iter+1] = theta.next
  theta.curr = theta.next 
}
# ------ End Simulation ------ # 

hist(theta.Samples[500:50000], breaks=seq(0,5,by=1))



# 3-2. select next sample : 후보 선택확률 보정 
theta.Samples = c(1:N)*0
theta.Samples[1] = theta.curr

# ------ Start Simulation ------ # 
prob.prop = c(.1, .1, .2, .3, .3)   # (1,2,3,4,5) 선택 될 확률 
for (iter in 1:(N-1)){
  theta.prop = sample(theta, 1, prob=prob.prop)     
  
  # 서로의 상대적 후보확률의 차이 보정 (새로운 이동확률 값 부여)
  alpha.star = prob.ratio(theta.prop, theta.curr, theta, prob) * 
    prob.ratio(theta.curr, theta.prop, theta, prob.prop)
  
  alpha = min(1, alpha.star)
  theta.next = ifelse(runif(1) < alpha, theta.prop, theta.curr)
  theta.Samples[iter+1] = theta.next
  theta.curr = theta.next 
}
# ------ End Simulation ------ # 

hist(theta.Samples[500:50000], breaks=seq(0,5,by=1))

