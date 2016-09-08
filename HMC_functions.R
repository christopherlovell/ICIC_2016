
## Hamiltonian Monte Carlo

H <- function(L,K){
  return(-log(L + K))
}

K <- function(u){
  (u %*% u) / 2
}




