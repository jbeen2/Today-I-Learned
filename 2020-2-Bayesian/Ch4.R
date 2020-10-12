prob.list <- rep(NA, 9)
for (i in 1:9) {
  prob.list[i] <- (factorial(20) / factorial(5) / factorial(15)) * (i*0.1)^5 * (1-i*0.1)^15
}

pi <- prob.list / sum(prob.list)
pi
0.0670+0.3667+ 0.3757+ 0.1568+ 0.0031 +0.0003


HPDgrid <- function(prob, level=0.95){
  prob.sort = sort(prob, decreasing=T)
  M = min(which(cumsum(prob.sort) >= level))
  height = prob.sort[M]
  HPD.index = which(prob>=height)
  HPD.level = sum(prob[HPD.index])
  res = list(index = HPD.index, level=HPD.level)
  return(res)
}


alpha = 0.05 ; level = 1-alpha
HPD = HPDgrid(pi, level)
HPDgrid.hat <- c(min(pi[HPD$index]), max(pi[HPD$index]))
HPDgrid.hat

plot(pi[HPD$index], rep(1, length(HPD$index)))
HPDgrid.hat[2] - HPDgrid.hat[1]
