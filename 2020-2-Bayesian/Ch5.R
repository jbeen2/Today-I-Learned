# Ch5 7p : 몬테카를로 추정, sample based estimation 
theta = rbeta(2000,16,26)
hist(theta, prob=T) # 상대도수 histogram
lines(density(theta) , lty=2) # 표본을 가지고 density estimate, 겹쳐 그리기 
theta.grid = seq(0.1, 0.7, length=100)
lines(theta.grid, dbeta(theta.grid, 16, 26), col=2) # true posterior density 
# 점선과 빨간색이 거의 비슷하다 -> 점선 사용해도 괜찮다! = 몬테카를로 (컴퓨터에서 난수 생성)
mean(theta)
var(theta)
# 정확도 더 높이고 싶으면 sample size(2000) 더 늘리기 

plot(density(theta), type="l")
CI = quantile(theta, prob=c(0.025, 0.975)) ; CI
abline(v=CI, col=3) # 평행한 직선 or 새로운 선 


# 16p
a=1 ; b=1 
theta = seq(0.01, 0.99, length=100)
plot(theta, dbeta(theta, 16, 26), type="l") # posterior 
lines(theta, dbeta(theta, 1, 1), col=2) 

a=2 ; b=10 # prior 변경 
lines(theta, dbeta(theta, 15+2, 25+10), col=3) 
lines(theta, dbeta(theta, 2, 10), col=4) # 상당히 많이 왼쪽으로 당겨짐 
15/40 # maximum 

lines(theta, dbeta(theta, 0.2, 1), col=5, lty=2) 
lines(theta, dbeta(theta, 15+0.2, 25+1), col=6) # 훨씬 덜 당겨진 것 
