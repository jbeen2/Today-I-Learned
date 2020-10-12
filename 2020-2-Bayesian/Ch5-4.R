# 1 p.133
a=1 ; b=1 
x=2 ; n=10 
theta = seq(0,1,length=1001)
ftheta = dbeta(theta, a+x, n-x+b)
prob = ftheta/sum(ftheta)

HPDgrid = function(prob, level=0.95) { 
  prob.sort = sort(prob, decreasing=T)
  M = min(which(cumsum(prob.sort) >= level))
  height = prob.sort[M]
  HPD.index = which(prob >= height)
  HPD.level = sum(prob[HPD.index])
  res = list(index=HPD.index, level=HPD.level)
  return(res)
}

HPD = HPDgrid(prob, 0.95)
HPD.grid = c(min(theta[HPD$index]), max(theta[HPD$index])) ; HPD.grid
