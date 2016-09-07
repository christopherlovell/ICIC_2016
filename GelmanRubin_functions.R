
B <- function(chains){
  
  theta_m <- unlist(lapply(chains, function(x) mean(x)))
  
  theta <- mean(theta_m)
  
  n <- length(unlist(chains[1]))
  M <- length(chains)
  
  return(n / (M-1) * sum(unlist(lapply(theta_m, function(x) (x - theta)^2))))
}

W <- function(chains){
  
  theta_m <- unlist(lapply(chains, function(x) mean(x)))
  
  return(mean(unlist(lapply(chains, function(x) 1 / (length(unlist(x)) - 1) * sum((unlist(x) - theta_m[1])^2)))))
}

Vhat <- function(chains){
  n <- length(unlist(chains[1]))
  M <- length(chains)
  
  return((n-1)/n * W(chains) + (M+1)/(n * M) * B(chains))
}